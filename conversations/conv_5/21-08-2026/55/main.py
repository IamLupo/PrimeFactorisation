#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# EXPERIMENT 616
# ==============================================================================
#
# CORRECTED RAW-RESIDUAL VALUATION LAW
#
# We work only with odd primes and odd semiprimes p*q.
#
# Frame A:
#
#     A = p - 3
#     B = q + 3
#
# Frame B:
#
#     A = p + 1
#     B = q - 3
#
# At level z=2:
#
#     a = A/2
#     b = B/2
#
# The normalized representation survives while a and b have equal parity.
#
# The experiment tests:
#
#     future_steps
#       =
#     min(v2(a), v2(b), v2(a-b)-1)
#
# and, after substituting the raw residuals:
#
#     depth
#       =
#     1 + min(v2(A), v2(B), v2(A-B)-1)
#
# A further simplification suggested by the valuation identity is:
#
#     depth = 1 + min(v2(A), v2(B)).
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
# PRIME SIEVE
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
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# SEMIPRIME GENERATION
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
    r = n % 4

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Expected odd semiprime n == 1 or 3 mod 4, "
        f"got n={n}, n mod 4={r}"
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
        x_num = q - p + 6
        y_num = p + q

    elif frame == "B":
        x_num = 3 * p - q + 6
        y_num = 3 * p + q

    else:
        raise ValueError(
            f"Unknown frame={frame}"
        )

    if x_num % 2 != 0 or y_num % 2 != 0:
        raise ArithmeticError(
            f"Non-integral X,Y: "
            f"frame={frame}, p={p}, q={q}, "
            f"numerators=({x_num},{y_num})"
        )

    return x_num // 2, y_num // 2


# ==============================================================================
# 2-ADIC VALUATION
# ==============================================================================

INF = 10**9


def v2(n: int) -> int:
    if n == 0:
        return INF

    n = abs(n)

    return (n & -n).bit_length() - 1


def fmt_v2_value(n: int) -> str:
    value = v2(n)
    return "inf" if value >= INF else str(value)


# ==============================================================================
# RAW LEVEL-2 RESIDUALS
# ==============================================================================

def raw_residuals(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":
        return (
            p - 3,
            q + 3,
        )

    if frame == "B":
        return (
            p + 1,
            q - 3,
        )

    raise ValueError(
        f"Unknown frame={frame}"
    )


# ==============================================================================
# NORMALIZED RESIDUALS AT LEVEL z
# ==============================================================================

def residuals_at_level(
    frame: str,
    p: int,
    q: int,
    z: int,
) -> tuple[int, int] | None:

    if z < 2:
        raise ValueError(
            f"Level must be >= 2, got z={z}"
        )

    modulus = 1 << (z - 1)

    A, B = raw_residuals(
        frame,
        p,
        q,
    )

    if A % modulus != 0:
        return None

    if B % modulus != 0:
        return None

    return (
        A // modulus,
        B // modulus,
    )


# ==============================================================================
# LOCAL X,Y FROM NORMALIZED RESIDUALS
# ==============================================================================

def local_xy_numerators(
    frame: str,
    a: int,
    b: int,
) -> tuple[int, int]:

    if frame == "A":
        return (
            b - a,
            a + b,
        )

    if frame == "B":
        return (
            3 * a - b,
            3 * a + b,
        )

    raise ValueError(
        f"Unknown frame={frame}"
    )


def local_xy_integral(
    frame: str,
    a: int,
    b: int,
) -> bool:

    x_num, y_num = local_xy_numerators(
        frame,
        a,
        b,
    )

    return (
        x_num % 2 == 0
        and
        y_num % 2 == 0
    )


# ==============================================================================
# OBSERVED FUTURE STEPS
# ==============================================================================

def observed_future_steps(
    frame: str,
    p: int,
    q: int,
    start_z: int,
) -> int:

    steps = 0

    for z in range(
        start_z,
        MAX_LEVEL,
    ):

        current = residuals_at_level(
            frame,
            p,
            q,
            z,
        )

        if current is None:
            break

        a, b = current

        if not local_xy_integral(
            frame,
            a,
            b,
        ):
            break

        nxt = residuals_at_level(
            frame,
            p,
            q,
            z + 1,
        )

        if nxt is None:
            break

        na, nb = nxt

        if not local_xy_integral(
            frame,
            na,
            nb,
        ):
            break

        steps += 1

    return steps


# ==============================================================================
# TEST 0
# ==============================================================================

def test_domain(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for state in states:

        if state.p % 2 == 0:
            failures += 1
            continue

        if state.q % 2 == 0:
            failures += 1
            continue

        if state.n % 2 == 0:
            failures += 1
            continue

        if state.n % 4 not in (1, 3):
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 1
# ==============================================================================

def test_corrected_raw_formula(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 1: CORRECTED RAW-RESIDUAL CLOSED FORM")
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

        predicted = (
            min(
                alpha,
                beta,
                gamma - 1,
            )
            - 1
        )

        actual = observed_future_steps(
            frame,
            state.p,
            state.q,
            2,
        )

        checked += 1

        if predicted != actual:

            failures += 1

            if failures <= EXAMPLE_LIMIT:
                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"A={A} B={B} "
                    f"v2A={fmt_v2_value(A)} "
                    f"v2B={fmt_v2_value(B)} "
                    f"v2diff={fmt_v2_value(A-B)} "
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

def test_simple_min_formula(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 2: SIMPLE RAW VALUATION DEPTH FORMULA")
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

        predicted_depth = (
            1
            + min(
                v2(A),
                v2(B),
            )
        )

        actual_depth = (
            2
            + observed_future_steps(
                frame,
                state.p,
                state.q,
                2,
            )
        )

        checked += 1

        if predicted_depth != actual_depth:

            failures += 1

            if failures <= EXAMPLE_LIMIT:
                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"A={A} B={B} "
                    f"v2A={fmt_v2_value(A)} "
                    f"v2B={fmt_v2_value(B)} "
                    f"actual={actual_depth} "
                    f"predicted={predicted_depth}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 3
# ==============================================================================

def test_global_depth_bridge(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 3: FACTOR DEPTH == GLOBAL DEPTH")
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

        global_depth = (
            1
            + min(
                v2(X),
                v2(Y),
            )
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        factor_depth = (
            1
            + min(
                v2(A),
                v2(B),
            )
        )

        checked += 1

        if global_depth != factor_depth:

            failures += 1

            if failures <= EXAMPLE_LIMIT:
                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"X={X} Y={Y} "
                    f"A={A} B={B} "
                    f"global={global_depth} "
                    f"factor={factor_depth}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 4
# ==============================================================================

def test_equal_unequal_v2(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 4: EQUAL VS UNEQUAL RAW v2")
    print("=" * 90)

    equal = 0
    unequal = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        if v2(A) == v2(B):
            equal += 1
        else:
            unequal += 1

        predicted = (
            1
            + min(
                v2(A),
                v2(B),
            )
        )

        actual = (
            2
            + observed_future_steps(
                frame,
                state.p,
                state.q,
                2,
            )
        )

        if predicted != actual:
            failures += 1

    print(
        f"equal valuation states={equal}"
    )

    print(
        f"unequal valuation states={unequal}"
    )

    print(
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 5
# ==============================================================================

def test_difference_identity(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 5: v2(A-B) IDENTITY")
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

        if alpha != beta:

            # For unequal 2-adic valuations:
            #
            #     v2(A-B) = min(v2(A),v2(B)).

            expected = min(
                alpha,
                beta,
            )

            if gamma != expected:

                failures += 1

                if failures <= EXAMPLE_LIMIT:
                    print(
                        f"    mismatch unequal "
                        f"n={state.n} "
                        f"A={A} B={B} "
                        f"v2A={fmt_v2_value(A)} "
                        f"v2B={fmt_v2_value(B)} "
                        f"v2diff={fmt_v2_value(A-B)}"
                    )

        else:

            # For equal valuations:
            #
            #     v2(A-B) > v2(A)
            #
            # unless A=B, where v2(A-B)=inf.

            if gamma <= alpha:

                failures += 1

                if failures <= EXAMPLE_LIMIT:
                    print(
                        f"    mismatch equal "
                        f"n={state.n} "
                        f"A={A} B={B} "
                        f"v2A=v2B={fmt_v2_value(A)} "
                        f"v2diff={fmt_v2_value(A-B)}"
                    )

        checked += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 6
# ==============================================================================

def test_all_levels(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 6: CURRENT-LEVEL CLOSED FORM")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        for z in range(
            2,
            MAX_LEVEL,
        ):

            residuals = residuals_at_level(
                frame,
                state.p,
                state.q,
                z,
            )

            if residuals is None:
                continue

            a, b = residuals

            if not local_xy_integral(
                frame,
                a,
                b,
            ):
                continue

            observed = observed_future_steps(
                frame,
                state.p,
                state.q,
                z,
            )

            predicted = min(
                v2(a),
                v2(b),
                v2(a - b) - 1,
            )

            checked += 1

            if observed != predicted:

                failures += 1

                if failures <= EXAMPLE_LIMIT:
                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={frame} "
                        f"a={a} b={b} "
                        f"actual={observed} "
                        f"predicted={predicted}"
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
    print("TEST 7: RESIDUAL HALVING RECURSION")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        for z in range(
            2,
            MAX_LEVEL - 1,
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

            if not local_xy_integral(
                frame,
                a,
                b,
            ):
                continue

            checked += 1

            if aa != a // 2 or bb != b // 2:

                failures += 1

                if failures <= EXAMPLE_LIMIT:
                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"({a},{b}) "
                        f"-> ({aa},{bb})"
                    )

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
    failures = 0

    for state in states:

        frame = frame_from_n(state.n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        predicted = (
            1
            + min(
                v2(A),
                v2(B),
            )
        )

        actual = (
            2
            + observed_future_steps(
                frame,
                state.p,
                state.q,
                2,
            )
        )

        distribution[actual] += 1

        if predicted != actual:
            failures += 1

    for depth in sorted(distribution):
        print(
            f"    depth={depth}: "
            f"{distribution[depth]}"
        )

    print()
    print(
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
        9,
        15,
        21,
        33,
        39,
        51,
        57,
        69,
        87,
        93,
        111,
        123,
        129,
        141,
        159,
        177,
        183,
        201,
        213,
        219,
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

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        raw_depth = (
            1
            + min(
                v2(A),
                v2(B),
            )
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
            f"    v2(X)={fmt_v2_value(X)} "
            f"v2(Y)={fmt_v2_value(Y)}"
        )

        print(
            f"    raw residuals="
            f"({A},{B})"
        )

        print(
            f"    v2(A)={fmt_v2_value(A)} "
            f"v2(B)={fmt_v2_value(B)} "
            f"v2(A-B)={fmt_v2_value(A-B)}"
        )

        print(
            f"    predicted depth={raw_depth}"
        )

        for z in range(
            2,
            min(raw_depth + 2, MAX_LEVEL),
        ):

            residuals = residuals_at_level(
                frame,
                state.p,
                state.q,
                z,
            )

            if residuals is None:
                break

            a, b = residuals

            x_num, y_num = local_xy_numerators(
                frame,
                a,
                b,
            )

            alive = local_xy_integral(
                frame,
                a,
                b,
            )

            print(
                f"    z={z:<2} "
                f"a={a:<8} "
                f"b={b:<8} "
                f"x={x_num}/2 "
                f"y={y_num}/2 "
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
Let the raw level-2 residuals be:

    FRAME A:

        A = p - 3
        B = q + 3

    FRAME B:

        A = p + 1
        B = q - 3.

At level z=2:

    a = A/2
    b = B/2.

The local factor representation is integral iff:

    a == b (mod 2).

For a represented normalized state:

    future_steps
      =
    min(
        v2(a),
        v2(b),
        v2(a-b)-1
    ).

Since:

    v2(a) = v2(A)-1
    v2(b) = v2(B)-1
    v2(a-b) = v2(A-B)-1,

we obtain:

    future_steps
      =
    min(
        v2(A),
        v2(B),
        v2(A-B)-1
    ) - 1.

The actual depth is therefore:

    depth
      =
    2 + future_steps

      =
    1 +
    min(
        v2(A),
        v2(B),
        v2(A-B)-1
    ).

Now use the standard 2-adic identity:

If:

    v2(A) != v2(B),

then:

    v2(A-B) = min(v2(A),v2(B)).

If:

    v2(A) = v2(B),

then:

    v2(A-B) > v2(A)

or A-B=0.

Therefore in BOTH cases:

    depth
      =
    1 + min(
        v2(A),
        v2(B)
    ).

Hence the candidate closed form is:

FRAME A:

    depth
      =
    1 + min(
        v2(p-3),
        v2(q+3)
    ).

FRAME B:

    depth
      =
    1 + min(
        v2(p+1),
        v2(q-3)
    ).

This experiment tests that formula directly against the
complete level-by-level representation.

The global coordinates must satisfy the identical result:

    depth
      =
    1 + min(
        v2(X),
        v2(Y)
    ).

Thus the desired bridge is:

    1 + min(
        v2(X),
        v2(Y)
    )

    =

    1 + min(
        v2(A),
        v2(B)
    ).

The level hierarchy would then reduce to the valuations of
two simple linear forms in the prime factors.
"""
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 616 START")
    print("=" * 90)
    print()
    print("CORRECTED RAW-RESIDUAL VALUATION LAW")
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

    test_domain(states)

    test_corrected_raw_formula(states)

    test_simple_min_formula(states)

    test_global_depth_bridge(states)

    test_equal_unequal_v2(states)

    test_difference_identity(states)

    test_all_levels(states)

    test_residual_halving(states)

    test_depth_distribution(states)

    print_examples(states)

    symbolic_summary()

    print("=" * 90)
    print("EXPERIMENT 616 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()