#!/usr/bin/env python3

import random
import sympy as sp


# ============================================================================
# EXPERIMENT 79
# STATIC r1/r2 -> n=p*q -> K,E,T,k,l,a,b
# ============================================================================

SEED = 1511464998

CASES = 1000

P_MIN = 10_000
P_MAX = 100_000

# Static prime values.
R1_VALUES = (
    29, 31, 37, 41, 43, 47, 53, 59, 61,
)

R2_VALUES = (
    67, 71, 73, 79, 83, 89, 97, 101, 103,
)


# ============================================================================
# PRIME GENERATION
# ============================================================================

def random_prime(rng, minimum, maximum):
    while True:
        x = rng.randint(minimum, maximum)
        if sp.isprime(x):
            return x


# ============================================================================
# PARAMETER CALCULATION
# ============================================================================

def calculate_parameters(p, q, r1, r2):
    """
    Decompose:

        p = r1*k + a
        q = r2*l + b

    with:

        0 <= a < r1
        0 <= b < r2

    Then define:

        R = r1*r2
        K = k*l
        T = floor(n/R)
        E = T-K
    """

    n = p * q
    R = r1 * r2

    # Euclidean decomposition of p and q.
    k, a = divmod(p, r1)
    l, b = divmod(q, r2)

    # Main quantities.
    K = k * l
    T = n // R
    E = T - K

    # ------------------------------------------------------------------------
    # Fundamental decomposition.
    # ------------------------------------------------------------------------

    assert p == r1 * k + a
    assert q == r2 * l + b

    assert 0 <= a < r1
    assert 0 <= b < r2

    # ------------------------------------------------------------------------
    # Expanded product.
    #
    # n = (r1*k+a)(r2*l+b)
    #
    #   = R*K
    #     + r1*k*b
    #     + r2*l*a
    #     + a*b
    # ------------------------------------------------------------------------

    cross = (
        r1 * k * b
        + r2 * l * a
        + a * b
    )

    assert n == R * K + cross

    # ------------------------------------------------------------------------
    # E is the quotient contribution of the cross term.
    # ------------------------------------------------------------------------

    expected_E = cross // R

    assert E == expected_E
    assert T == K + E

    # ------------------------------------------------------------------------
    # The remaining sub-cell residue is < R.
    # ------------------------------------------------------------------------

    remainder = cross - R * E

    assert 0 <= remainder < R

    # ------------------------------------------------------------------------
    # Equivalent quotient decomposition.
    #
    # n = R*T + remainder
    # ------------------------------------------------------------------------

    assert n == R * T + remainder

    # ------------------------------------------------------------------------
    # Direct definition of T.
    # ------------------------------------------------------------------------

    assert T == n // R
    assert n % R == remainder

    return (
        n,
        R,
        T,
        K,
        E,
        k,
        l,
        a,
        b,
    )


# ============================================================================
# MAIN
# ============================================================================

def main():
    rng = random.Random(SEED)

    print(
        "n	p	q	r1	r2	R	T	K	E	k	l	a	b"
    )

    for case_index in range(1, CASES + 1):

        # --------------------------------------------------------------------
        # Generate two different prime factors.
        # --------------------------------------------------------------------

        while True:
            p = random_prime(rng, P_MIN, P_MAX)
            q = random_prime(rng, P_MIN, P_MAX)

            if p != q:
                break

        # --------------------------------------------------------------------
        # Select static r1/r2.
        # --------------------------------------------------------------------

        r1 = rng.choice(R1_VALUES)
        r2 = rng.choice(R2_VALUES)

        (
            n,
            R,
            T,
            K,
            E,
            k,
            l,
            a,
            b,
        ) = calculate_parameters(
            p,
            q,
            r1,
            r2,
        )

        print(
            f"{n}	"
            f"{p}	"
            f"{q}	"
            f"{r1}	"
            f"{r2}	"
            f"{R}	"
            f"{T}	"
            f"{K}	"
            f"{E}	"
            f"{k}	"
            f"{l}	"
            f"{a}	"
            f"{b}"
        )


if __name__ == "__main__":
    main()