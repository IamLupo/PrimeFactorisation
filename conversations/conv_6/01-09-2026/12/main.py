#!/usr/bin/env python3

import math
import random
import time


# ==========================================================================================
# EXPERIMENT 83
# MULTISCALE HIDDEN-STATE TRAJECTORY
#
# For one fixed:
#
#     n = p*q
#
# we generate:
#
#     v_0 = n
#     v_1 = floor(n/2)
#     v_2 = floor(n/4)
#     ...
#
# BUT:
#
#     v is NEVER factored.
#
# v is used ONLY to select:
#
#     r1 = previous_prime(sqrt(v))
#     r2 = next_prime(sqrt(v))
#
# Then all hidden parameters are computed against the ORIGINAL n,p,q:
#
#     p = r1*k + a
#     q = r2*l + b
#
#     K = k*l
#
#     R = r1*r2
#
#     Q = floor(n/R)
#
#     E = Q-K
#
#     c1 = floor(k*b/r2)
#     c2 = floor(l*a/r1)
#
#     beta  = (k*b) % r2
#     alpha = (l*a) % r1
#
#     c3 = floor((r1*beta + r2*alpha + a*b)/R)
#
#     E_check = c1+c2+c3
#
# We print everything so transitions between levels can be studied.
# ==========================================================================================


SEED = 1511464998
#random.seed(SEED)


# ------------------------------------------------------------------------------------------
# 64-bit primality
# ------------------------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    # Write n-1 = d * 2^s
    d = n - 1
    s = 0

    while (d & 1) == 0:
        d >>= 1
        s += 1

    # Deterministic for unsigned 64-bit integers.
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

        witness = True

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                witness = False
                break

        if witness:
            return False

    return True


def previous_prime(x: int) -> int:
    x = int(x)

    if x <= 2:
        return 2

    if x == 3:
        return 2

    if x % 2 == 0:
        x -= 1

    while x >= 3:
        if is_prime(x):
            return x
        x -= 2

    return 2


def next_prime(x: int) -> int:
    x = int(x)

    if x <= 2:
        return 2

    if x % 2 == 0:
        x += 1

    while True:
        if is_prime(x):
            return x
        x += 2


# ------------------------------------------------------------------------------------------
# Hidden-state calculation
# ------------------------------------------------------------------------------------------

def calculate_hidden_state(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
):
    R = r1 * r2

    # Quotient coordinates.
    k = p // r1
    l = q // r2

    # Residue coordinates.
    a = p % r1
    b = q % r2

    # Product quotient coordinate.
    K = k * l

    # n/R quotient.
    Q = n // R

    # Error/carry coordinate.
    E = Q - K

    # Carry decomposition.
    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // R

    E_carry = c1 + c2 + c3

    # Exact reconstruction identity.
    reconstructed_n = (
        R * K
        + r1 * k * b
        + r2 * l * a
        + a * b
    )

    return {
        "R": R,
        "Q": Q,

        "K": K,
        "E": E,

        "k": k,
        "l": l,

        "a": a,
        "b": b,

        "c1": c1,
        "c2": c2,
        "c3": c3,

        "alpha": alpha,
        "beta": beta,

        "E_carry": E_carry,
        "reconstructed_n": reconstructed_n,

        "a_norm_num": a,
        "a_norm_den": r1,
        "b_norm_num": b,
        "b_norm_den": r2,

        "p_recovered": r1 * k + a,
        "q_recovered": r2 * l + b,
    }


# ------------------------------------------------------------------------------------------
# Transition metrics
# ------------------------------------------------------------------------------------------

def add_transition_fields(prev, cur):
    if prev is None:
        return

    integer_fields = (
        "R",
        "Q",
        "K",
        "E",
        "k",
        "l",
        "a",
        "b",
        "c1",
        "c2",
        "c3",
        "alpha",
        "beta",
    )

    for field in integer_fields:
        cur["d_" + field] = cur[field] - prev[field]

    # Normalized residue coordinates.
    cur["a_ratio"] = cur["a"] / cur["r1"]
    cur["b_ratio"] = cur["b"] / cur["r2"]

    prev_a_ratio = prev["a"] / prev["r1"]
    prev_b_ratio = prev["b"] / prev["r2"]

    cur["d_a_ratio"] = cur["a_ratio"] - prev_a_ratio
    cur["d_b_ratio"] = cur["b_ratio"] - prev_b_ratio


# ------------------------------------------------------------------------------------------
# Generate semiprime
# ------------------------------------------------------------------------------------------

def generate_prime_near(target: int) -> int:
    """
    Generate a prime near target.
    """
    if target < 10:
        target = 10

    while True:
        candidate = random.randint(1, target)

        if candidate < 3:
            continue

        if candidate % 2 == 0:
            candidate += 1

        if is_prime(candidate):
            return candidate


def generate_semiprime(scale: int):
    """
    Generate two roughly sqrt(scale)-sized primes so that:

        n = p*q

    is near the requested scale.
    """

    root = int(math.isqrt(scale))

    p = generate_prime_near(root)
    q = generate_prime_near(root)

    if(p > q):
        p, q = q, p

    return p * q, p, q


# ------------------------------------------------------------------------------------------
# Print level
# ------------------------------------------------------------------------------------------

def print_level(
    level,
    n,
    v,
    sqrt_v,
    r1,
    r2,
    state,
):
    print(f"LEVEL {level}")
    print(f"    v                   = {v:,}")
    print(f"    sqrt(v)             = {sqrt_v:,}")

    print(f"    r1                  = {r1:,}")
    print(f"    r2                  = {r2:,}")
    print(f"    R=r1*r2             = {state['R']:,}")
    print(f"    n/R                 = {n / state['R']:.15f}")
    print(f"    floor(n/R)          = {state['Q']:,}")

    print()
    print(f"    K                   = {state['K']:,}")
    print(f"    E                   = {state['E']:,}")
    print(f"    k                   = {state['k']:,}")
    print(f"    l                   = {state['l']:,}")
    print(f"    a                   = {state['a']:,}")
    print(f"    b                   = {state['b']:,}")

    print()
    print(f"    c1                  = {state['c1']:,}")
    print(f"    c2                  = {state['c2']:,}")
    print(f"    c3                  = {state['c3']:,}")
    print(f"    c1+c2+c3            = {state['E_carry']:,}")

    print()
    print(f"    alpha=(l*a)%r1     = {state['alpha']:,}")
    print(f"    beta =(k*b)%r2     = {state['beta']:,}")

    print()
    print(f"    a/r1                = {state['a'] / r1:.12f}")
    print(f"    b/r2                = {state['b'] / r2:.12f}")

    print()
    print(f"    r1*k+a              = {state['p_recovered']:,}")
    print(f"    r2*l+b              = {state['q_recovered']:,}")
    print(f"    p*q                 = {state['p_recovered'] * state['q_recovered']:,}")
    print(f"    matches n           = {state['reconstructed_n'] == n}")

    if "d_K" in state:
        print()
        print("    TRANSITION FROM PREVIOUS LEVEL")

        print(f"        dK              = {state['d_K']:,}")
        print(f"        dE              = {state['d_E']:,}")
        print(f"        dk              = {state['d_k']:,}")
        print(f"        dl              = {state['d_l']:,}")
        print(f"        da              = {state['d_a']:,}")
        print(f"        db              = {state['d_b']:,}")

        print(f"        dc1             = {state['d_c1']:,}")
        print(f"        dc2             = {state['d_c2']:,}")
        print(f"        dc3             = {state['d_c3']:,}")

        print(f"        d(alpha)        = {state['d_alpha']:,}")
        print(f"        d(beta)         = {state['d_beta']:,}")

        print(f"        d(a/r1)        = {state['d_a_ratio']:.12f}")
        print(f"        d(b/r2)        = {state['d_b_ratio']:.12f}")

    print()
    print("-" * 92)


# ------------------------------------------------------------------------------------------
# Analyze one n
# ------------------------------------------------------------------------------------------

def run_number(
    case_index: int,
    scale_label: str,
    scale: int,
    levels: int,
):
    n, p, q = generate_semiprime(scale)

    print()
    print("=" * 92)
    print(f"CASE {case_index}")
    print("=" * 92)

    print(f"scale target         = {scale_label}")
    print(f"n                    = {n:,}")
    print(f"true p               = {p:,}")
    print(f"true q               = {q:,}")
    print(f"p*q == n             = {p*q == n}")
    print()

    previous = None
    level = 0

    #for level in range(levels):
    while n >> level >= 1:
        # IMPORTANT:
        # v changes, but p and q do NOT.
        v = n >> level

        if v <= 4:
            break

        sqrt_v = math.isqrt(v)

        #r1 = previous_prime(sqrt_v)
        #r2 = next_prime(sqrt_v)

        r1 = sqrt_v
        r2 = sqrt_v
        

        state = calculate_hidden_state(
            n=n,
            p=p,
            q=q,
            r1=r1,
            r2=r2,
        )

        # Attach r1/r2 because normalized and transition analysis
        # uses them.
        state["r1"] = r1
        state["r2"] = r2

        add_transition_fields(previous, state)

        print_level(
            level=level,
            n=n,
            v=v,
            sqrt_v=sqrt_v,
            r1=r1,
            r2=r2,
            state=state,
        )

        previous = state
        level += 1


# ------------------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------------------

def main():
    print("=" * 92)
    print("START EXPERIMENT 83")
    print("MULTISCALE HIDDEN-STATE TRAJECTORY")
    print("=" * 92)

    print()
    print("configuration")
    print("    levels / n         = 20")
    print("    cases / scale      = 3")
    print("    scales             = ['1e9', '1e12', '1e16']")
    print("    v                  = floor(n / 2^level)")
    print("    v is NEVER factored")
    print("    r1/r2 from sqrt(v)")
    print("    K,E,k,l,a,b from ORIGINAL n,p,q")
    print(f"    seed               = {SEED}")
    print()

    cases_per_scale = 3
    levels = 20

    scales = [
        ("1e9", 10**9),
        ("1e12", 10**12),
        ("1e16", 10**16),
    ]

    start = time.perf_counter()

    case_index = 1

    for scale_label, scale in scales:
        print()
        print("=" * 92)
        print(f"SCALE {scale_label}")
        print("=" * 92)

        for _ in range(cases_per_scale):
            run_number(
                case_index=case_index,
                scale_label=scale_label,
                scale=scale,
                levels=levels,
            )
            case_index += 1

    elapsed = time.perf_counter() - start

    print()
    print("=" * 92)
    print("FINISHED EXPERIMENT 83")
    print("=" * 92)
    print(f"total runtime = {elapsed:.6f}s")
    print()


if __name__ == "__main__":
    main()