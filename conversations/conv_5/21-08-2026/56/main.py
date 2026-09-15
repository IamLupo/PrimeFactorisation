#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from fractions import Fraction


# ==============================================================================
# EXPERIMENT 617
# ==============================================================================
#
# EXACT PIECEWISE 2-ADIC DEPTH LAW
#
# Experiment 616 showed:
#
#   simple depth = 1 + min(v2(A), v2(B))
#
# is WRONG exactly when the two valuations differ.
#
# The observed law from Experiment 611 is:
#
#   alpha = v2(A)
#   beta  = v2(B)
#
#   depth =
#
#       min(alpha,beta) + 1    if alpha == beta
#
#       min(alpha,beta)        if alpha != beta
#
# Equivalently:
#
#   depth = min(alpha,beta) + [alpha == beta]
#
# where [condition] is 1 when true and 0 otherwise.
#
# This experiment:
#
#   1. Verifies the piecewise raw-residual depth law.
#   2. Verifies the equivalent v2(A-B) formulation.
#   3. Verifies the global X,Y depth law.
#   4. Proves the factor/global valuation correspondence.
#   5. Searches for a gcd-based formulation.
#   6. Tests the law level-by-level.
#   7. Produces a compact classification table.
#
# Domain:
#
#   odd primes only
#   p <= q
#   odd semiprimes n=pq
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250
MAX_LEVEL = 24
EXAMPLE_LIMIT = 30


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


# ==============================================================================
# SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:
    if limit < 3:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            sieve[start : limit + 1 : p] = b"\x00" * count

    return [
        p for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# SEMIPRIMES
# ==============================================================================

def generate_states(primes: list[int]) -> list[State]:
    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append(
                State(
                    n=p * q,
                    p=p,
                    q=q,
                )
            )

    return states


# ==============================================================================
# FRAME
# ==============================================================================

def frame_from_n(n: int) -> str:
    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Invalid odd semiprime residue: "
        f"n={n}, n mod 4={r}"
    )


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":
        xn = q - p + 6
        yn = p + q

    elif frame == "B":
        xn = 3 * p - q + 6
        yn = 3 * p + q

    else:
        raise ValueError(frame)

    if xn % 2 or yn % 2:
        raise ArithmeticError(
            f"Non-integral global coordinates: "
            f"frame={frame}, p={p}, q={q}"
        )

    return xn // 2, yn // 2


# ==============================================================================
# RAW FACTOR RESIDUALS
# ==============================================================================

def raw_residuals(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":
        return p - 3, q + 3

    if frame == "B":
        return p + 1, q - 3

    raise ValueError(frame)


# ==============================================================================
# 2-ADIC VALUATION
# ==============================================================================

INF = 10**9


def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def vstr(x: int) -> str:
    y = v2(x)
    return "inf" if y >= INF else str(y)


# ==============================================================================
# PIECEWISE DEPTH LAW
# ==============================================================================

def predicted_factor_depth(
    frame: str,
    p: int,
    q: int,
) -> int:

    A, B = raw_residuals(
        frame,
        p,
        q,
    )

    alpha = v2(A)
    beta = v2(B)

    return min(alpha, beta) + int(alpha == beta)


# ==============================================================================
# GLOBAL DEPTH
# ==============================================================================

def predicted_global_depth(
    X: int,
    Y: int,
) -> int:

    return 1 + min(
        v2(X),
        v2(Y),
    )


# ==============================================================================
# RESIDUALS AT LEVEL
# ==============================================================================

def residuals_at_level(
    frame: str,
    p: int,
    q: int,
    z: int,
) -> tuple[int, int] | None:

    modulus = 1 << (z - 1)

    A, B = raw_residuals(
        frame,
        p,
        q,
    )

    if A % modulus:
        return None

    if B % modulus:
        return None

    return A // modulus, B // modulus


# ==============================================================================
# LOCAL REPRESENTATION
# ==============================================================================

def local_xy_numerators(
    frame: str,
    a: int,
    b: int,
) -> tuple[int, int]:

    if frame == "A":
        return b - a, a + b

    if frame == "B":
        return 3 * a - b, 3 * a + b

    raise ValueError(frame)


def local_integral(
    frame: str,
    a: int,
    b: int,
) -> bool:

    xn, yn = local_xy_numerators(
        frame,
        a,
        b,
    )

    return xn % 2 == 0 and yn % 2 == 0


# ==============================================================================
# OBSERVED DEPTH
# ==============================================================================

def observed_depth(
    frame: str,
    p: int,
    q: int,
) -> int:

    depth = 1

    for z in range(2, MAX_LEVEL + 1):

        r = residuals_at_level(
            frame,
            p,
            q,
            z,
        )

        if r is None:
            break

        a, b = r

        if not local_integral(
            frame,
            a,
            b,
        ):
            break

        depth = z

    return depth


# ==============================================================================
# TEST 1
# ==============================================================================

def test_piecewise_depth(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 1: EXACT PIECEWISE RAW-RESIDUAL DEPTH LAW")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        actual = observed_depth(
            frame,
            state.p,
            state.q,
        )

        predicted = predicted_factor_depth(
            frame,
            state.p,
            state.q,
        )

        checked += 1

        if actual != predicted:

            failures += 1

            if failures <= EXAMPLE_LIMIT:

                A, B = raw_residuals(
                    frame,
                    state.p,
                    state.q,
                )

                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"A={A} B={B} "
                    f"v2A={vstr(A)} "
                    f"v2B={vstr(B)} "
                    f"actual={actual} "
                    f"predicted={predicted}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 2
# ==============================================================================

def test_v2_difference_formula(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 2: EQUIVALENT v2(A-B) FORMULA")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)
        gamma = v2(A - B)

        # Exact formula:
        #
        # depth = 1 + min(
        #     alpha,
        #     beta,
        #     gamma - 1
        # )
        #
        # This is valid for the represented base state.

        predicted = (
            1
            + min(
                alpha,
                beta,
                gamma - 1,
            )
        )

        actual = observed_depth(
            frame,
            state.p,
            state.q,
        )

        checked += 1

        if predicted != actual:

            failures += 1

            if failures <= EXAMPLE_LIMIT:
                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"alpha={vstr(A)} "
                    f"beta={vstr(B)} "
                    f"gamma={vstr(A-B)} "
                    f"actual={actual} "
                    f"predicted={predicted}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 3
# ==============================================================================

def test_global_depth(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 3: RAW FACTOR DEPTH == GLOBAL X,Y DEPTH")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        factor_depth = predicted_factor_depth(
            frame,
            state.p,
            state.q,
        )

        global_depth = predicted_global_depth(
            X,
            Y,
        )

        checked += 1

        if factor_depth != global_depth:

            failures += 1

            if failures <= EXAMPLE_LIMIT:

                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"X={X} Y={Y} "
                    f"factor={factor_depth} "
                    f"global={global_depth}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 4
# ==============================================================================

def test_piecewise_classes(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 4: EQUAL / UNEQUAL v2 CLASSIFICATION")
    print("=" * 90)

    table = Counter()

    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)

        key = (
            frame,
            "equal" if alpha == beta else "unequal",
            min(alpha, beta),
        )

        table[key] += 1

        predicted = (
            min(alpha, beta)
            + int(alpha == beta)
        )

        actual = observed_depth(
            frame,
            state.p,
            state.q,
        )

        if predicted != actual:
            failures += 1

    for key in sorted(table, key=str):
        print(
            f"    {key} -> "
            f"{table[key]}"
        )

    print()
    print(
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 5
# ==============================================================================

def test_valuation_identity(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 5: 2-ADIC DIFFERENCE IDENTITY")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)
        gamma = v2(A - B)

        checked += 1

        if alpha != beta:

            expected = min(
                alpha,
                beta,
            )

            if gamma != expected:
                failures += 1

        else:

            if gamma <= alpha:
                failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 6
# ==============================================================================

def test_level_formula(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 6: EXACT CURRENT-LEVEL DEPTH FORMULA")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        for z in range(
            2,
            MAX_LEVEL,
        ):

            r = residuals_at_level(
                frame,
                state.p,
                state.q,
                z,
            )

            if r is None:
                continue

            a, b = r

            if not local_integral(
                frame,
                a,
                b,
            ):
                continue

            alpha = v2(a)
            beta = v2(b)

            predicted_future = (
                min(alpha, beta)
                + int(alpha == beta)
            )

            # Because current level z is represented, the number
            # of *additional* represented levels is one less:
            #
            #   future = depth - z
            #
            # and the next failure occurs after consuming
            # the current valuation layer.

            predicted_depth = (
                z
                + predicted_future
                - 1
            )

            actual_depth = observed_depth(
                frame,
                state.p,
                state.q,
            )

            checked += 1

            if predicted_depth != actual_depth:

                failures += 1

                if failures <= EXAMPLE_LIMIT:
                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"a={a} b={b} "
                        f"actual={actual_depth} "
                        f"predicted={predicted_depth}"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 7
# ==============================================================================

def test_residual_halving(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 7: RESIDUAL HALVING")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        for z in range(
            2,
            MAX_LEVEL,
        ):

            r0 = residuals_at_level(
                frame,
                state.p,
                state.q,
                z,
            )

            r1 = residuals_at_level(
                frame,
                state.p,
                state.q,
                z + 1,
            )

            if r0 is None or r1 is None:
                continue

            a, b = r0
            aa, bb = r1

            if not local_integral(
                frame,
                a,
                b,
            ):
                continue

            checked += 1

            if aa != a // 2 or bb != b // 2:
                failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 8
# ==============================================================================

def test_depth_distribution(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 8: DEPTH DISTRIBUTION")
    print("=" * 90)

    distribution = Counter()

    for state in states:

        frame = frame_from_n(state.n)

        depth = observed_depth(
            frame,
            state.p,
            state.q,
        )

        distribution[depth] += 1

    for depth in sorted(distribution):
        print(
            f"    depth={depth}: "
            f"{distribution[depth]}"
        )

    print()


# ==============================================================================
# TEST 9
# ==============================================================================

def test_gcd_formulation(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 9: GCD / 2-ADIC FORMULATION")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        g = abs(__import__("math").gcd(A, B))

        # v2(gcd(A,B)) = min(v2(A),v2(B))
        #
        # The piecewise law therefore becomes:
        #
        #   depth = v2(gcd(A,B)) + 1
        #             if v2(A)=v2(B)
        #
        #   depth = v2(gcd(A,B))
        #             otherwise.

        gval = v2(g)

        equal = (
            v2(A) == v2(B)
        )

        predicted = (
            gval + int(equal)
        )

        actual = observed_depth(
            frame,
            state.p,
            state.q,
        )

        checked += 1

        if predicted != actual:
            failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states: list[State],
) -> None:

    requested = [
        9, 15, 21, 33, 39,
        51, 57, 69, 87, 93,
        111, 123, 129, 141, 159,
        177, 183, 201, 213, 219,
    ]

    lookup = {
        state.n: state
        for state in states
    }

    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    for n in requested:

        state = lookup.get(n)

        if state is None:
            continue

        frame = frame_from_n(state.n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)

        depth = observed_depth(
            frame,
            state.p,
            state.q,
        )

        print()
        print(
            f"n={n} "
            f"p={state.p} "
            f"q={state.q}"
        )

        print(
            f"    frame={frame}"
        )

        print(
            f"    X={X} Y={Y}"
        )

        print(
            f"    raw residuals="
            f"({A},{B})"
        )

        print(
            f"    v2(A)={vstr(A)} "
            f"v2(B)={vstr(B)} "
            f"equal={alpha == beta}"
        )

        print(
            f"    "
            f"min_v2={min(alpha,beta)} "
            f"piecewise_depth="
            f"{min(alpha,beta)+int(alpha==beta)} "
            f"observed_depth={depth}"
        )

        for z in range(
            2,
            min(depth + 2, MAX_LEVEL),
        ):

            r = residuals_at_level(
                frame,
                state.p,
                state.q,
                z,
            )

            if r is None:
                break

            a, b = r

            alive = local_integral(
                frame,
                a,
                b,
            )

            print(
                f"    z={z:<2} "
                f"a={a:<8} "
                f"b={b:<8} "
                f"bits=({a & 1},{b & 1}) "
                f"alive={alive}"
            )

            if not alive:
                break

    print()


# ==============================================================================
# SYMBOLIC SUMMARY
# ==============================================================================

def symbolic_summary() -> None:

    print("=" * 90)
    print("SYMBOLIC SUMMARY")
    print("=" * 90)

    print(
r"""
Let:

    A,B = raw factor residuals at level 2.

FRAME A:

    A = p - 3
    B = q + 3

FRAME B:

    A = p + 1
    B = q - 3.

Define:

    alpha = v2(A)
    beta  = v2(B).

The experimentally exact depth law is:

    depth =
        min(alpha,beta) + 1,   alpha = beta

        min(alpha,beta),       alpha != beta.

Equivalently:

    depth =
        min(alpha,beta)
        + [alpha = beta].

The reason is the standard 2-adic identity:

    alpha != beta
        =>
    v2(A-B)=min(alpha,beta)

whereas:

    alpha = beta
        =>
    v2(A-B)>alpha
    or A=B.

At a represented normalized level:

    a = A / 2^(z-1)
    b = B / 2^(z-1)

and the representation condition is:

    a == b (mod 2).

The next residual pair is:

    (a,b) -> (a/2,b/2)

and the process terminates exactly when the two
valuations stop matching.

GLOBAL COORDINATES:

    depth =
        1 + min(v2(X),v2(Y)).

Thus the key identity is:

    1 + min(v2(X),v2(Y))
      =
    min(v2(A),v2(B))
      + [v2(A)=v2(B)].

This is the correct replacement for the false
formula:

    1 + min(v2(A),v2(B)).

The extra +1 appears only in the equal-valuation branch.

For the two historical frames:

    A-frame:
        A = p-3
        B = q+3

    B-frame:
        A = p+1
        B = q-3.

Therefore the factor-side hierarchy is governed by the
2-adic valuations of two linear forms in the factors,
with a single equality flag deciding whether the terminal
level receives one additional step.
"""
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 617 START")
    print("=" * 90)
    print()
    print("EXACT PIECEWISE 2-ADIC DEPTH LAW")
    print()

    print("[1] PRIME SIEVE")

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )
    print()

    print("[2] SEMIPRIME GENERATION")

    states = generate_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    print("[3] GLOBAL STATES")

    print(
        f"    global states={len(states)}"
    )
    print()

    test_piecewise_depth(states)

    test_v2_difference_formula(states)

    test_global_depth(states)

    test_piecewise_classes(states)

    test_valuation_identity(states)

    test_level_formula(states)

    test_residual_halving(states)

    test_depth_distribution(states)

    test_gcd_formulation(states)

    print_examples(states)

    symbolic_summary()

    print("=" * 90)
    print("EXPERIMENT 617 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
