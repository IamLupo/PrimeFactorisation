#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 583
# ==============================================================================
# NORMALIZED COORDINATE / PARENT -> CHILD RESIDUE TRANSITION
#
# NO FILES
# NO WEB
#
# We use:
#
#     u = y - x
#     v = y + x
#
# and the user's level-dependent formulas:
#
#     k = 2^(z-1)
#
# CASE A:
#
#     n mod 2^z = 2^z - 1
#
#     p = k*y - k*x + 3
#     q = k*y + k*x - 3
#
# therefore:
#
#     p = k*u + 3
#     q = k*v - 3
#
#
# CASE B:
#
#     otherwise
#
#     p = (k*y + k*x - 3)/3
#     q = k*y - k*x + 3
#
# therefore:
#
#     3p = k*v - 3
#     q  = k*u + 3
#
#
# Main hypothesis
# ---------------
#
# Because:
#
#     k = 2^(z-1)
#
# and:
#
#     k^2 is divisible by 2^z,
#
# the next residue should depend on very little information about u,v.
#
# We explicitly test:
#
#     (u mod 2, v mod 2)
#
# and determine whether that pair predicts the child residue exactly.
#
# We also derive the symbolic residue transition for CASE A and CASE B.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# CONFIG
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 8

SHOW_EXAMPLES = True
MAX_EXAMPLES = 12


# ==============================================================================
# SAMPLE
# ==============================================================================

@dataclass(frozen=True)
class Sample:

    n: int
    p: int
    q: int

    z: int
    k: int
    modulus: int

    branch: str

    x: int
    y: int

    u: int
    v: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int):

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = isqrt(limit)

    for p in range(
        2,
        root + 1,
    ):

        if not sieve[p]:
            continue

        start = p * p

        count = (
            (limit - start) // p
        ) + 1

        sieve[
            start:limit + 1:p
        ] = b"\x00" * count

    return sieve


# ==============================================================================
# SEMIPRIMES
# ==============================================================================

def generate_semiprimes(
    max_n,
    sieve,
):

    primes = [
        p
        for p in range(
            3,
            max_n + 1,
            2,
        )
        if sieve[p]
    ]

    result = []

    for i, p in enumerate(primes):

        if p * p > max_n:
            break

        max_q = max_n // p

        for q in primes[i:]:

            if q > max_q:
                break

            result.append(
                (
                    p * q,
                    p,
                    q,
                )
            )

    return result


# ==============================================================================
# SOLVE CASE A
# ==============================================================================

def solve_case_a(
    p,
    q,
    k,
):

    # p = ku + 3
    # q = kv - 3

    if (p - 3) % k != 0:
        return None

    if (q + 3) % k != 0:
        return None

    u = (
        p - 3
    ) // k

    v = (
        q + 3
    ) // k

    # x = (v-u)/2
    # y = (v+u)/2

    if (v - u) % 2 != 0:
        return None

    if (v + u) % 2 != 0:
        return None

    x = (
        v - u
    ) // 2

    y = (
        v + u
    ) // 2

    return x, y, u, v


# ==============================================================================
# SOLVE CASE B
# ==============================================================================

def solve_case_b(
    p,
    q,
    k,
):

    # 3p = kv - 3
    # q  = ku + 3

    if (3 * p + 3) % k != 0:
        return None

    if (q - 3) % k != 0:
        return None

    v = (
        3 * p + 3
    ) // k

    u = (
        q - 3
    ) // k

    if (v - u) % 2 != 0:
        return None

    if (v + u) % 2 != 0:
        return None

    x = (
        v - u
    ) // 2

    y = (
        v + u
    ) // 2

    return x, y, u, v


# ==============================================================================
# GENERATE LEVEL SAMPLE
# ==============================================================================

def make_level_sample(
    n,
    p,
    q,
    z,
):

    k = 1 << (
        z - 1
    )

    modulus = 1 << z

    residue = n % modulus

    if residue == modulus - 1:

        solved = solve_case_a(
            p,
            q,
            k,
        )

        if solved is None:
            return None

        x, y, u, v = solved

        branch = "A"

        check_p = (
            k * u
            + 3
        )

        check_q = (
            k * v
            - 3
        )

    else:

        solved = solve_case_b(
            p,
            q,
            k,
        )

        if solved is None:
            return None

        x, y, u, v = solved

        branch = "B"

        numerator = (
            k * v
            - 3
        )

        if numerator % 3 != 0:
            return None

        check_p = numerator // 3

        check_q = (
            k * u
            + 3
        )

    if (
        check_p != p
        or check_q != q
    ):
        raise RuntimeError(
            "Factor reconstruction failed"
        )

    if (
        p * q != n
    ):
        raise RuntimeError(
            "n != p*q"
        )

    return Sample(
        n=n,
        p=p,
        q=q,
        z=z,
        k=k,
        modulus=modulus,
        branch=branch,
        x=x,
        y=y,
        u=u,
        v=v,
    )


# ==============================================================================
# BUILD DATA
# ==============================================================================

def generate_data(
    semiprimes,
):

    data = defaultdict(list)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        for n, p, q in semiprimes:

            sample = make_level_sample(
                n,
                p,
                q,
                z,
            )

            if sample is not None:

                data[z].append(
                    sample
                )

    return data


# ==============================================================================
# SYMBOLIC CASE-A RESIDUE
# ==============================================================================

def symbolic_case_a(
    z,
):

    k = 1 << (
        z - 1
    )

    modulus = 1 << z

    # n = (ku+3)(kv-3)
    #
    #   = k^2 uv - 3ku + 3kv - 9
    #
    # modulo 2^z:
    #
    # k^2uv vanishes.
    #
    # 3k(v-u) has only the top bit of k surviving.

    print()
    print(
        f"CASE A SYMBOLIC z={z}"
    )

    print(
        f"    k = {k}"
    )

    print(
        f"    modulus = {modulus}"
    )

    print()
    print(
        "    n = (ku+3)(kv-3)"
    )

    print(
        "      = k^2uv + 3k(v-u) - 9"
    )

    print(
        "    mod 2^z:"
    )

    print(
        "      k^2uv -> 0"
    )

    print(
        "      3k(v-u) depends only on"
    )

    print(
        "      (v-u) mod 2"
    )

    if z >= 2:

        residue_even = (
            -9
        ) % modulus

        residue_odd = (
            -9
            + 3 * k
        ) % modulus

        print()
        print(
            "    if (v-u) mod 2 = 0:"
        )

        print(
            f"        n = {residue_even} mod {modulus}"
        )

        print(
            "    if (v-u) mod 2 = 1:"
        )

        print(
            f"        n = {residue_odd} mod {modulus}"
        )


# ==============================================================================
# SYMBOLIC CASE-B RESIDUE
# ==============================================================================

def symbolic_case_b(
    z,
):

    k = 1 << (
        z - 1
    )

    modulus = 1 << z

    print()
    print(
        f"CASE B SYMBOLIC z={z}"
    )

    print(
        f"    k = {k}"
    )

    print(
        f"    modulus = {modulus}"
    )

    # n = (ku+3)(kv-3)/3
    #
    #   = [k^2uv + 3k(v-u) - 9]/3
    #
    # The modular behavior requires the numerator to be divisible by 3.
    #
    # Since 2^z is invertible mod 3, the reduction can be examined through
    # the parity state together with the divisibility condition.

    print()
    print(
        "    n = (ku+3)(kv-3)/3"
    )

    print(
        "      = [k^2uv + 3k(v-u) - 9]/3"
    )

    print()
    print(
        "    The k^2uv term vanishes modulo 2^z"
    )

    print(
        "    after accounting for the division by 3,"
    )

    print(
        "    because 3 is invertible modulo 2^z."
    )

    inverse_3 = pow(
        3,
        -1,
        modulus,
    )

    even_residue = (
        (
            -9
        )
        * inverse_3
    ) % modulus

    odd_residue = (
        (
            -9
            + 3 * k
        )
        * inverse_3
    ) % modulus

    print()
    print(
        f"    inverse(3) mod {modulus} = "
        f"{inverse_3}"
    )

    print()
    print(
        "    if (v-u) mod 2 = 0:"
    )

    print(
        f"        n = {even_residue} mod {modulus}"
    )

    print(
        "    if (v-u) mod 2 = 1:"
    )

    print(
        f"        n = {odd_residue} mod {modulus}"
    )


# ==============================================================================
# EMPIRICAL PARITY TEST
# ==============================================================================

def parity_rule_test(
    samples,
):

    print()
    print("=" * 90)
    print(
        "EMPIRICAL (u-v) PARITY TEST"
    )
    print("=" * 90)

    for branch in (
        "A",
        "B",
    ):

        print()
        print(
            f"BRANCH {branch}"
        )

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            subset = [
                s
                for s in samples[z]
                if s.branch == branch
            ]

            if not subset:
                print(
                    f"    z={z}: no samples"
                )
                continue

            table = defaultdict(
                set
            )

            modulus = 1 << z

            for s in subset:

                parity = (
                    s.v - s.u
                ) & 1

                child = (
                    s.n
                    % (
                        1 << (z + 1)
                    )
                )

                table[
                    parity
                ].add(
                    child
                )

            deterministic = all(
                len(values) == 1
                for values
                in table.values()
            )

            print(
                f"    z={z:<2} "
                f"parent mod={modulus:<4} "
                f"deterministic={deterministic}"
            )

            for parity in sorted(
                table
            ):

                print(
                    f"        parity={parity} "
                    f"child residues="
                    f"{sorted(table[parity])}"
                )


# ==============================================================================
# DIRECT CHILD BIT TEST
# ==============================================================================

def child_bit(
    n,
    z,
):

    return (
        n >> z
    ) & 1


def child_bit_parity_test(
    samples,
):

    print()
    print("=" * 90)
    print(
        "NEXT BIT FROM (v-u) PARITY"
    )
    print("=" * 90)

    for branch in (
        "A",
        "B",
    ):

        print()
        print(
            f"BRANCH {branch}"
        )

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            subset = [
                s
                for s in samples[z]
                if s.branch == branch
            ]

            if not subset:
                continue

            mapping = defaultdict(
                set
            )

            for s in subset:

                parity = (
                    s.v - s.u
                ) & 1

                bit = child_bit(
                    s.n,
                    z,
                )

                mapping[
                    parity
                ].add(
                    bit
                )

            exact = all(
                len(values) == 1
                for values
                in mapping.values()
            )

            print(
                f"    z={z:<2} "
                f"exact={exact}"
            )

            for parity in sorted(
                mapping
            ):

                print(
                    f"        parity={parity} "
                    f"bits={sorted(mapping[parity])}"
                )


# ==============================================================================
# FULL (u mod 2, v mod 2) TEST
# ==============================================================================

def uv_parity_test(
    samples,
):

    print()
    print("=" * 90)
    print(
        "FULL (u mod 2, v mod 2) STATE TEST"
    )
    print("=" * 90)

    for branch in (
        "A",
        "B",
    ):

        print()
        print(
            f"BRANCH {branch}"
        )

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            subset = [
                s
                for s in samples[z]
                if s.branch == branch
            ]

            if not subset:
                continue

            table = defaultdict(
                set
            )

            for s in subset:

                state = (
                    s.u & 1,
                    s.v & 1,
                )

                bit = child_bit(
                    s.n,
                    z,
                )

                table[
                    state
                ].add(
                    bit
                )

            exact = all(
                len(values) == 1
                for values in table.values()
            )

            print(
                f"    z={z:<2} "
                f"exact={exact}"
            )

            for state in sorted(
                table
            ):

                print(
                    f"        "
                    f"(u,v) parity={state} "
                    f"bits={sorted(table[state])}"
                )


# ==============================================================================
# DOMAIN SHRINKING IN NORMALIZED COORDINATES
# ==============================================================================

def normalized_domain(
    samples,
):

    print()
    print("=" * 90)
    print(
        "NORMALIZED DOMAIN"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print()
        print(
            f"z={z}"
        )

        for branch in (
            "A",
            "B",
        ):

            subset = [
                s
                for s in samples[z]
                if s.branch == branch
            ]

            if not subset:
                continue

            us = {
                s.u
                for s in subset
            }

            vs = {
                s.v
                for s in subset
            }

            uv = {
                (
                    s.u,
                    s.v,
                )
                for s in subset
            }

            print(
                f"    branch {branch}: "
                f"states={len(subset):<7} "
                f"u={len(us):<7} "
                f"v={len(vs):<7} "
                f"uv={len(uv):<7} "
                f"u_range=[{min(us)},"
                f"{max(us)}] "
                f"v_range=[{min(vs)},"
                f"{max(vs)}]"
            )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    samples,
):

    if not SHOW_EXAMPLES:
        return

    print()
    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        subset = samples[z]

        print()
        print(
            f"z={z}"
        )

        for s in subset[
            :MAX_EXAMPLES
        ]:

            print(
                f"    n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"x={s.x} "
                f"y={s.y} "
                f"u={s.u} "
                f"v={s.v} "
                f"branch={s.branch}"
            )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 583 START"
    )
    print("=" * 90)

    print()
    print(
        "NORMALIZED COORDINATE TRANSITION"
    )

    # --------------------------------------------------------------------------
    # Generate factor states.
    # --------------------------------------------------------------------------

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = prime_sieve(
        MAX_N
    )

    print(
        "[2] SEMIPRIME GENERATION"
    )

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        f"    semiprimes={len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Build level data.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILDING LEVEL STATES"
    )

    data = generate_data(
        semiprimes
    )

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print(
            f"    z={z} "
            f"states={len(data[z])}"
        )

    # --------------------------------------------------------------------------
    # Symbolic derivation.
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print(
        "SYMBOLIC TRANSITIONS"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        symbolic_case_a(
            z
        )

        symbolic_case_b(
            z
        )

    # --------------------------------------------------------------------------
    # Empirical tests.
    # --------------------------------------------------------------------------

    child_bit_parity_test(
        data
    )

    uv_parity_test(
        data
    )

    parity_rule_test(
        data
    )

    # --------------------------------------------------------------------------
    # Normalized domains.
    # --------------------------------------------------------------------------

    normalized_domain(
        data
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    print_examples(
        data
    )

    # --------------------------------------------------------------------------
    # Final.
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print(
        "FINAL AUDIT"
    )
    print("=" * 90)

    print(
        """
The experiment tests whether the modular hierarchy collapses to
a tiny normalized state.

Definitions:

    u = y-x
    v = y+x

Case A:

    p = ku+3
    q = kv-3

Case B:

    3p = kv-3
    q  = ku+3

with:

    k = 2^(z-1).

Because k^2 is divisible by 2^z, the product n loses the
full uv dependence modulo 2^z.

The critical test is whether:

    (u mod 2, v mod 2)

or even:

    (v-u) mod 2

determines the next branch bit exactly.

If that holds repeatedly, the modular hierarchy has a very small
state transition despite p and q becoming much larger.

That is precisely the kind of state reduction that could be useful
for a SAT encoding.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 583 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
