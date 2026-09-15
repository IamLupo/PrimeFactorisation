#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 124
# Carry-state elimination via cross-row / cross-column rank-1 recovery
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 124")
print("Carry-state elimination via cross-row / cross-column recovery")
print("=" * 72)


# ------------------------------------------------------------------------
# GRID
# ------------------------------------------------------------------------

R1 = [17, 43]
R2 = [19, 47]

R00 = R1[0] * R2[0]
R01 = R1[0] * R2[1]
R10 = R1[1] * R2[0]
R11 = R1[1] * R2[1]

print()
print("R1 =", R1)
print("R2 =", R2)
print()
print("R00 =", R00)
print("R01 =", R01)
print("R10 =", R10)
print("R11 =", R11)


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    # Enough for the sizes used in this experiment.
    bases = (
        2, 325, 9375, 28178,
        450775, 9780504, 1795265022
    )

    for a in bases:
        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n

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
        x |= (1 << (bits - 1))
        x |= 1

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# SEMIPRIME GENERATOR
# ------------------------------------------------------------------------

def random_semiprime(bits):
    pb = bits // 2
    qb = bits - pb

    while True:
        p = random_prime(pb)
        q = random_prime(qb)

        if p == q:
            continue

        n = p * q

        if n.bit_length() == bits:
            return n, min(p, q), max(p, q)


# ------------------------------------------------------------------------
# CARRY CALCULATION
# ------------------------------------------------------------------------

def carries(p, q, r1, r2):
    k, a = divmod(p, r1)
    ell, b = divmod(q, r2)

    c1 = (k * b) // r2
    c2 = (ell * a) // r1

    beta = (k * b) % r2
    alpha = (ell * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    E = c1 + c2 + c3

    return {
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "E": E,
    }


# ------------------------------------------------------------------------
# QUOTIENT
# ------------------------------------------------------------------------

def quotient(n, r1, r2):
    return n // (r1 * r2)


# ------------------------------------------------------------------------
# EXACT CARRY GRID
# ------------------------------------------------------------------------

def build_true_grid(n, p, q):
    cells = {}

    for i, r1 in enumerate(R1):
        for j, r2 in enumerate(R2):
            C = carries(p, q, r1, r2)

            Q = quotient(n, r1, r2)

            cells[(i, j)] = {
                "Q": Q,
                "E": C["E"],
                "K": Q - C["E"],
                "k": C["k"],
                "ell": C["ell"],
                "a": C["a"],
                "b": C["b"],
            }

    return cells


# ------------------------------------------------------------------------
# DIRECT SQRT SCAN
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
# NEW METHOD
#
# For a given row pair:
#
#   K00 = k0*l0
#   K10 = k1*l0
#
# and:
#
#   p = r10*k1 + a1
#     = r00*k0 + a0
#
# Therefore:
#
#   r00*K00 - r10*K10
#       = l0 * (a1-a0)
#
# giving
#
#   l0 =
#       [r00*K00 - r10*K10] / (a1-a0)
#
# Once l0 is known:
#
#   k0 = K00/l0
#   p  = r00*k0+a0
#
# Then q follows from n/p.
#
# The difficulty is that K = Q-E and E is unknown.
#
# We therefore enumerate bounded carry totals rather than p.
# ------------------------------------------------------------------------


def recovery_from_carries(
    n,
    e00,
    e10,
    a0,
    a1,
    b0,
):
    r00 = R1[0]
    r10 = R1[1]

    Q00 = quotient(n, R1[0], R2[0])
    Q10 = quotient(n, R1[1], R2[0])

    K00 = Q00 - e00
    K10 = Q10 - e10

    if K00 <= 0 or K10 <= 0:
        return None

    d = a1 - a0

    # ------------------------------------------------------------
    # Same residue.
    #
    # If a1 == a0 then the equation reduces to:
    #
    #     r00*K00 = r10*K10
    #
    # and l0 cannot be obtained from this row pair.
    # ------------------------------------------------------------

    if d == 0:
        return None

    numerator = r00 * K00 - r10 * K10

    if numerator <= 0:
        return None

    if numerator % d != 0:
        return None

    ell0 = numerator // d

    if ell0 <= 0:
        return None

    if K00 % ell0 != 0:
        return None

    k0 = K00 // ell0

    p = r00 * k0 + a0

    if p <= 1:
        return None

    if p >= n:
        return None

    if n % p != 0:
        return None

    q = n // p

    # ------------------------------------------------------------
    # Check q residue for the chosen b0.
    # ------------------------------------------------------------

    if q % R2[0] != b0:
        return None

    return p, q, k0, ell0


# ------------------------------------------------------------------------
# CARRY BOUNDS
#
# Since
#
#   c1 < k
#   c2 < ell
#   c3 <= 2
#
# we have
#
#   E <= k + ell + 2.
#
# For the experiment we derive a conservative bound using
# sqrt(n).
# ------------------------------------------------------------------------

def global_carry_bound(n):
    s = math.isqrt(n)

    # Extremely conservative.
    #
    # k < p/r1 <= sqrt(n)/min(r1)
    # ell < q/r2 <= sqrt(n)/min(r2)
    #
    # Use max possible values.
    max_k = s // min(R1) + 1
    max_l = s // min(R2) + 1

    return max_k + max_l + 2


# ------------------------------------------------------------------------
# NEW SEARCH
# ------------------------------------------------------------------------

def carry_state_search(n):
    Q00 = quotient(n, R1[0], R2[0])
    Q10 = quotient(n, R1[1], R2[0])

    E_MAX = global_carry_bound(n)

    print()
    print("Q00 =", Q00)
    print("Q10 =", Q10)
    print("global E bound =", E_MAX)

    tests = 0
    residue_states = 0
    divisibility_hits = 0
    recovered = []

    # ------------------------------------------------------------
    # Enumerate residue state:
    #
    #   a0 in [0,r10)
    #   a1 in [0,r11)
    #   b0 in [0,r20)
    #
    # The second q residue is not needed for the direct
    # reconstruction. It is used later as verification.
    # ------------------------------------------------------------

    for a0 in range(R1[0]):
        for a1 in range(R1[1]):

            if a0 == a1:
                continue

            for b0 in range(R2[0]):

                residue_states += 1

                # ------------------------------------------------
                # Instead of blindly scanning E from zero to E_MAX,
                # derive a narrow relation between e00 and e10.
                #
                # We have:
                #
                # numerator =
                #   r00(Q00-e00) - r10(Q10-e10)
                #
                # and this must equal:
                #
                #   ell0*(a1-a0).
                #
                # Thus it must be divisible by d.
                #
                # We still scan e00/e10 here, but only over
                # residue-compatible values. This experiment tests
                # whether the resulting state space is manageable.
                # ------------------------------------------------

                d = a1 - a0

                # Keep the smaller carry dimension on the outside.
                for e00 in range(E_MAX + 1):

                    # Quick upper/lower constraints on numerator.
                    base = (
                        R1[0] * (Q00 - e00)
                        - R1[1] * Q10
                    )

                    if base <= 0:
                        continue

                    # Need:
                    #
                    # base + R1[1]*e10
                    #
                    # to be positive and divisible by d.
                    #
                    # Search e10 values, but use modular stepping.
                    #
                    # d can be negative.
                    if d > 0:
                        first = (-base) % d
                    else:
                        dd = -d
                        first = (-base) % dd

                    step = abs(d)

                    # If d == +-1, every e10 is possible.
                    if step == 0:
                        continue

                    e10 = first

                    while e10 <= E_MAX:

                        tests += 1

                        numerator = (
                            R1[0] * (Q00 - e00)
                            - R1[1] * (Q10 - e10)
                        )

                        if numerator > 0:
                            if numerator % d == 0:
                                divisibility_hits += 1

                                result = recovery_from_carries(
                                    n,
                                    e00,
                                    e10,
                                    a0,
                                    a1,
                                    b0,
                                )

                                if result is not None:
                                    p, q, k0, ell0 = result

                                    candidate = (
                                        p,
                                        q,
                                        a0,
                                        a1,
                                        b0,
                                        e00,
                                        e10,
                                        k0,
                                        ell0,
                                    )

                                    if candidate not in recovered:
                                        recovered.append(candidate)

                        e10 += step

    return {
        "E_MAX": E_MAX,
        "residue_states": residue_states,
        "carry_tests": tests,
        "divisibility_hits": divisibility_hits,
        "recovered": recovered,
    }


# ------------------------------------------------------------------------
# RUN EXPERIMENT
# ------------------------------------------------------------------------

random.seed(124)

TEST_BITS = [
    30,
    35,
    41,
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
    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("TRUE GAP")
    print("q-p =", true_q - true_p)

    true_grid = build_true_grid(n, true_p, true_q)

    print()
    print("TRUE CARRY GRID")

    for key in sorted(true_grid):
        cell = true_grid[key]

        print(
            f"cell {key}: "
            f"Q={cell['Q']} "
            f"E={cell['E']} "
            f"K={cell['K']} "
            f"k={cell['k']} "
            f"ell={cell['ell']} "
            f"a={cell['a']} "
            f"b={cell['b']}"
        )

    # ------------------------------------------------------------
    # Verify true state is inside E bound.
    # ------------------------------------------------------------

    E_MAX = global_carry_bound(n)

    max_true_E = max(
        cell["E"]
        for cell in true_grid.values()
    )

    print()
    print("MAX TRUE E =", max_true_E)
    print("E bound    =", E_MAX)

    if max_true_E <= E_MAX:
        print("TRUE E STATE INSIDE BOUND = True")
    else:
        print("TRUE E STATE INSIDE BOUND = False")

    # ------------------------------------------------------------
    # Direct control
    # ------------------------------------------------------------

    print()
    print("DIRECT SQRT(n) SCAN")

    t0 = time.perf_counter()

    p_direct, q_direct, tested_direct = direct_factor(n)

    t1 = time.perf_counter()

    direct_time = t1 - t0

    direct_ok = (
        p_direct * q_direct == n
        and {
            p_direct,
            q_direct
        } == {
            true_p,
            true_q
        }
    )

    print("p tested =", tested_direct)
    print("runtime  =", f"{direct_time:.6f}", "s")
    print("correct  =", direct_ok)

    # ------------------------------------------------------------
    # New carry-state method
    # ------------------------------------------------------------

    print()
    print("CARRY-STATE ELIMINATION")

    t0 = time.perf_counter()

    result = carry_state_search(n)

    t1 = time.perf_counter()

    carry_time = t1 - t0

    recovered = result["recovered"]

    correct_recovery = any(
        {candidate[0], candidate[1]} == {
            true_p,
            true_q
        }
        for candidate in recovered
    )

    print()
    print("E_MAX =", result["E_MAX"])
    print("residue states =", result["residue_states"])
    print("carry tests    =", result["carry_tests"])
    print("divisibility hits =", result["divisibility_hits"])
    print("recovered candidates =", len(recovered))
    print("runtime =", f"{carry_time:.6f}", "s")
    print("TRUE FACTORIZATION RECOVERED =", correct_recovery)

    # ------------------------------------------------------------
    # Print recoveries
    # ------------------------------------------------------------

    if recovered:
        print()
        print("RECOVERED STATES")

        for candidate in recovered[:20]:
            (
                p,
                q,
                a0,
                a1,
                b0,
                e00,
                e10,
                k0,
                ell0,
            ) = candidate

            exact = (
                {p, q}
                == {true_p, true_q}
            )

            print(
                "p =", p,
                "q =", q,
                "a0 =", a0,
                "a1 =", a1,
                "b0 =", b0,
                "E00 =", e00,
                "E10 =", e10,
                "k0 =", k0,
                "ell0 =", ell0,
                "EXACT =", exact
            )

    # ------------------------------------------------------------
    # Timing comparison
    # ------------------------------------------------------------

    print()
    print("TIMING")

    if carry_time > 0:
        print(
            "direct / carry-state =",
            f"{direct_time / carry_time:.3f}x"
        )

    # ------------------------------------------------------------
    # Final structural verification
    # ------------------------------------------------------------

    if correct_recovery:

        exact_candidate = next(
            candidate
            for candidate in recovered
            if {candidate[0], candidate[1]}
            == {true_p, true_q}
        )

        ep, eq = exact_candidate[:2]

        print()
        print("VERIFYING FULL CARRY GRID")

        grid_ok = True

        for i, r1 in enumerate(R1):
            for j, r2 in enumerate(R2):

                actual = carries(ep, eq, r1, r2)
                expected = true_grid[(i, j)]

                if (
                    actual["E"] != expected["E"]
                    or actual["k"] != expected["k"]
                    or actual["ell"] != expected["ell"]
                    or actual["a"] != expected["a"]
                    or actual["b"] != expected["b"]
                ):
                    grid_ok = False

        print("carry grid verified =", grid_ok)

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 124
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 124")
print("=" * 72)
