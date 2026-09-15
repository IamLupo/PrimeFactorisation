#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 614
==========================================================================================

EXACT 2-ADIC RESIDUAL LADDER — SAFE INTEGER/DEAD STATE HANDLING

Previous failure:

    ArithmeticError:
        Non-integral local coordinates: frame=B a=1 b=0

This is NOT an error.

For frame B:

    x = (3a-b)/2
    y = (3a+b)/2

so:

    a=1, b=0

gives:

    x=3/2
    y=3/2.

That simply means the factor residual pair does not represent an
integer normalized A/B state at that level.

This experiment therefore distinguishes:

    REPRESENTABLE
    DEAD / NON-INTEGRAL
    INVALID

instead of raising exceptions.

The core tests are:

    1. odd-semiprime domain
    2. residual representation criterion
    3. exact local-coordinate integrality
    4. one-step survival
    5. multi-level residual ladder
    6. valuation prediction
    7. residual/global X,Y crosswalk
    8. global depth agreement
    9. residual halving recursion
   10. minimal residual modulus

==========================================================================================
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isqrt
from fractions import Fraction


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250
MAX_LEVEL = 24
MAX_EXAMPLES = 20


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
            sieve[start : limit + 1 : p] = (
                b"\x00" * count
            )

    # IMPORTANT:
    # only odd primes are admitted.
    return [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_odd_semiprimes(
    primes: list[int],
) -> list[State]:

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
        f"Invalid odd-semiprime frame state: "
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
        raise ValueError(
            f"Unknown frame={frame}"
        )

    if (xn & 1) or (yn & 1):
        raise ArithmeticError(
            f"Global X,Y not integral: "
            f"frame={frame} p={p} q={q}"
        )

    return xn // 2, yn // 2


# ==============================================================================
# RESIDUALS
# ==============================================================================

def residuals_at_level(
    frame: str,
    p: int,
    q: int,
    z: int,
) -> tuple[int, int] | None:

    M = 1 << (z - 1)

    if frame == "A":

        rp = p - 3
        rq = q + 3

    elif frame == "B":

        rp = p + 1
        rq = q - 3

    else:
        raise ValueError(
            f"Unknown frame={frame}"
        )

    if rp % M != 0:
        return None

    if rq % M != 0:
        return None

    return (
        rp // M,
        rq // M,
    )


# ==============================================================================
# LOCAL COORDINATES
# ==============================================================================

def local_xy_fractional(
    frame: str,
    a: int,
    b: int,
) -> tuple[Fraction, Fraction]:

    if frame == "A":

        x = Fraction(b - a, 2)
        y = Fraction(a + b, 2)

    elif frame == "B":

        x = Fraction(3 * a - b, 2)
        y = Fraction(3 * a + b, 2)

    else:
        raise ValueError(
            f"Unknown frame={frame}"
        )

    return x, y


def local_xy_integer(
    frame: str,
    a: int,
    b: int,
) -> tuple[int, int] | None:

    x, y = local_xy_fractional(
        frame,
        a,
        b,
    )

    if x.denominator != 1:
        return None

    if y.denominator != 1:
        return None

    return (
        x.numerator,
        y.numerator,
    )


# ==============================================================================
# REPRESENTABILITY
# ==============================================================================

def residual_representation(
    frame: str,
    a: int,
    b: int,
) -> bool:

    return (
        local_xy_integer(
            frame,
            a,
            b,
        )
        is not None
    )


# ==============================================================================
# ONE-STEP SURVIVAL
# ==============================================================================

def residual_one_step(
    frame: str,
    a: int,
    b: int,
) -> bool:

    """
    Current state must first be representable.

    Then the next level must also be representable.

    Since the next residuals are a/2,b/2, this is equivalent to:

        a,b divisible by 2
        and
        (a/2,b/2) same parity.

    The exact test is performed directly rather than assumed.
    """

    if not residual_representation(
        frame,
        a,
        b,
    ):
        return False

    if a % 2 != 0:
        return False

    if b % 2 != 0:
        return False

    return residual_representation(
        frame,
        a // 2,
        b // 2,
    )


# ==============================================================================
# NEXT RESIDUALS
# ==============================================================================

def next_residuals(
    a: int,
    b: int,
) -> tuple[int, int] | None:

    if a % 2 != 0:
        return None

    if b % 2 != 0:
        return None

    return (
        a // 2,
        b // 2,
    )


# ==============================================================================
# 2-ADIC VALUATION
# ==============================================================================

INF_V2 = 10**9


def v2(n: int) -> int:

    if n == 0:
        return INF_V2

    n = abs(n)

    return (
        n & -n
    ).bit_length() - 1


def fmt_v2(n: int) -> str:

    value = v2(n)

    if value >= INF_V2:
        return "inf"

    return str(value)


# ==============================================================================
# GLOBAL LOCAL COORDINATES
# ==============================================================================

def global_local_xy(
    X: int,
    Y: int,
    z: int,
) -> tuple[int, int] | None:

    M = 1 << (z - 1)

    if X % M != 0:
        return None

    if Y % M != 0:
        return None

    return (
        X // M,
        Y // M,
    )


# ==============================================================================
# ACTUAL FUTURE DEPTH
# ==============================================================================

def actual_future_steps(
    frame: str,
    p: int,
    q: int,
    z: int,
) -> int:

    steps = 0
    level = z

    while level < MAX_LEVEL:

        current = residuals_at_level(
            frame,
            p,
            q,
            level,
        )

        if current is None:
            break

        a, b = current

        if not residual_representation(
            frame,
            a,
            b,
        ):
            break

        nxt = residuals_at_level(
            frame,
            p,
            q,
            level + 1,
        )

        if nxt is None:
            break

        na, nb = nxt

        if not residual_representation(
            frame,
            na,
            nb,
        ):
            break

        steps += 1
        level += 1

    return steps


# ==============================================================================
# CLOSED FORM FUTURE DEPTH
# ==============================================================================

def predicted_future_steps(
    frame: str,
    a: int,
    b: int,
) -> int:

    if not residual_representation(
        frame,
        a,
        b,
    ):
        return -1

    """
    A currently representable pair satisfies:

        a == b mod 2.

    Each additional level requires the next pair
    to remain representable.

    Let:

        va = v2(a)
        vb = v2(b)
        vd = v2(a-b)

    The maximum number of successful transitions is:

        min(
            va,
            vb,
            vd - 1
        )

    provided the current state is already representable.
    """

    va = v2(a)
    vb = v2(b)
    vd = v2(a - b)

    return max(
        0,
        min(
            va,
            vb,
            vd - 1,
        ),
    )


# ==============================================================================
# TEST 0 — DOMAIN
# ==============================================================================

def test_domain(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for state in states:

        if state.p == 2:
            failures += 1

        if state.q == 2:
            failures += 1

        if (state.n & 1) == 0:
            failures += 1

        if state.n % 4 not in (1, 3):
            failures += 1

    print(
        f"states={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()


# ==============================================================================
# TEST 1 — CURRENT REPRESENTATION
# ==============================================================================

def test_current_representation(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 1: CURRENT RESIDUAL REPRESENTATION")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        for z in range(
            2,
            MAX_LEVEL + 1,
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

            local = local_xy_integer(
                frame,
                a,
                b,
            )

            predicted = residual_representation(
                frame,
                a,
                b,
            )

            actual = (
                local is not None
            )

            checked += 1

            if predicted != actual:

                failures += 1

                if failures <= MAX_EXAMPLES:
                    xf, yf = local_xy_fractional(
                        frame,
                        a,
                        b,
                    )

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={frame} "
                        f"a={a} "
                        f"b={b} "
                        f"x={xf} "
                        f"y={yf}"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()


# ==============================================================================
# TEST 2 — ONE STEP
# ==============================================================================

def test_one_step(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 2: EXACT ONE-STEP RESIDUAL SURVIVAL")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        for z in range(
            2,
            MAX_LEVEL,
        ):

            current = residuals_at_level(
                frame,
                state.p,
                state.q,
                z,
            )

            if current is None:
                continue

            a, b = current

            if not residual_representation(
                frame,
                a,
                b,
            ):
                continue

            next_r = residuals_at_level(
                frame,
                state.p,
                state.q,
                z + 1,
            )

            actual = (
                next_r is not None
                and residual_representation(
                    frame,
                    *next_r,
                )
            )

            predicted = residual_one_step(
                frame,
                a,
                b,
            )

            checked += 1

            if actual != predicted:

                failures += 1

                if failures <= MAX_EXAMPLES:
                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={frame} "
                        f"a={a} "
                        f"b={b} "
                        f"actual={actual} "
                        f"predicted={predicted}"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()


# ==============================================================================
# TEST 3 — MULTI-LEVEL LADDER
# ==============================================================================

def test_multi_level_ladder(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 3: MULTI-LEVEL 2-ADIC LADDER")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        for z in range(
            2,
            MAX_LEVEL,
        ):

            current = residuals_at_level(
                frame,
                state.p,
                state.q,
                z,
            )

            if current is None:
                continue

            a, b = current

            if not residual_representation(
                frame,
                a,
                b,
            ):
                continue

            actual = actual_future_steps(
                frame,
                state.p,
                state.q,
                z,
            )

            predicted = predicted_future_steps(
                frame,
                a,
                b,
            )

            checked += 1

            if actual != predicted:

                failures += 1

                if failures <= MAX_EXAMPLES:
                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={frame} "
                        f"a={a} "
                        f"b={b} "
                        f"v2a={fmt_v2(a)} "
                        f"v2b={fmt_v2(b)} "
                        f"v2diff={fmt_v2(a-b)} "
                        f"actual={actual} "
                        f"predicted={predicted}"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()


# ==============================================================================
# TEST 4 — RESIDUAL HALVING
# ==============================================================================

def test_residual_halving(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 4: RESIDUAL HALVING RECURSION")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

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
            ap, bp = r1

            if not residual_one_step(
                frame,
                a,
                b,
            ):
                continue

            checked += 1

            if (
                ap != a // 2
                or bp != b // 2
            ):

                failures += 1

                if failures <= MAX_EXAMPLES:
                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"current=({a},{b}) "
                        f"expected=({a//2},{b//2}) "
                        f"actual=({ap},{bp})"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()


# ==============================================================================
# TEST 5 — GLOBAL X,Y CROSSWALK
# ==============================================================================

def test_global_crosswalk(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 5: GLOBAL X,Y <-> RESIDUAL COORDINATES")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

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

            local_residual = local_xy_integer(
                frame,
                a,
                b,
            )

            local_global = global_local_xy(
                X,
                Y,
                z,
            )

            if local_residual is None:
                # Dead state is allowed.
                if local_global is not None:
                    failures += 1

                    if failures <= MAX_EXAMPLES:
                        print(
                            f"    mismatch "
                            f"n={state.n} "
                            f"z={z} "
                            f"residual dead "
                            f"but global coordinates exist "
                            f"{local_global}"
                        )

                checked += 1
                continue

            checked += 1

            if local_residual != local_global:

                failures += 1

                if failures <= MAX_EXAMPLES:
                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={frame} "
                        f"a={a} "
                        f"b={b} "
                        f"residual_xy={local_residual} "
                        f"global_xy={local_global}"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()


# ==============================================================================
# TEST 6 — DEPTH
# ==============================================================================

def test_depth(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 6: EXACT GLOBAL DEPTH")
    print("=" * 90)

    failures = 0
    distribution = Counter()

    for state in states:

        frame = frame_from_n(
            state.n
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        predicted = (
            1
            + min(
                v2(X),
                v2(Y),
            )
        )

        actual = (
            2
            + actual_future_steps(
                frame,
                state.p,
                state.q,
                2,
            )
        )

        distribution[actual] += 1

        if predicted != actual:

            failures += 1

            if failures <= MAX_EXAMPLES:
                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"X={X} "
                    f"Y={Y} "
                    f"predicted={predicted} "
                    f"actual={actual}"
                )

    for depth in sorted(
        distribution
    ):
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
# TEST 7 — MINIMAL MODULUS
# ==============================================================================

def test_minimal_modulus(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 7: MINIMAL RESIDUAL MODULUS")
    print("=" * 90)

    for future in range(
        0,
        8,
    ):

        modulus = 1 << (future + 1)

        mapping: dict[
            tuple[str, int, int],
            set[bool],
        ] = {}

        for state in states:

            frame = frame_from_n(
                state.n
            )

            r = residuals_at_level(
                frame,
                state.p,
                state.q,
                2,
            )

            if r is None:
                continue

            a, b = r

            actual = (
                actual_future_steps(
                    frame,
                    state.p,
                    state.q,
                    2,
                )
                >= future
            )

            key = (
                frame,
                a % modulus,
                b % modulus,
            )

            mapping.setdefault(
                key,
                set(),
            ).add(actual)

        ambiguous = sum(
            1
            for outcomes in mapping.values()
            if len(outcomes) > 1
        )

        print(
            f"future={future} "
            f"modulus={modulus} "
            f"states={len(mapping)} "
            f"ambiguous={ambiguous}"
        )

    print()


# ==============================================================================
# TEST 8 — FACTOR RESIDUAL BIT STREAM
# ==============================================================================

def test_bit_stream(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 8: FACTOR RESIDUAL BIT STREAM")
    print("=" * 90)

    failures = 0
    checked = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

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

            if r0 is None:
                continue

            a, b = r0

            if not residual_representation(
                frame,
                a,
                b,
            ):
                continue

            r1 = residuals_at_level(
                frame,
                state.p,
                state.q,
                z + 1,
            )

            if r1 is None:
                continue

            ap, bp = r1

            checked += 1

            expected = (
                a // 2,
                b // 2,
            )

            if (ap, bp) != expected:

                failures += 1

                if failures <= MAX_EXAMPLES:
                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"r=({a},{b}) "
                        f"expected={expected} "
                        f"actual=({ap},{bp})"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()


# ==============================================================================
# EXAMPLES
# ==============================================================================

def examples(
    states: list[State],
) -> None:

    wanted = {
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
    }

    lookup = {
        state.n: state
        for state in states
        if state.n in wanted
    }

    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    for n in sorted(wanted):

        state = lookup.get(n)

        if state is None:
            continue

        frame = frame_from_n(
            state.n
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        print()
        print(
            f"n={state.n} "
            f"p={state.p} "
            f"q={state.q}"
        )

        print(
            f"    frame={frame} "
            f"X={X} Y={Y}"
        )

        print(
            f"    v2X={fmt_v2(X)} "
            f"v2Y={fmt_v2(Y)}"
        )

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
                break

            a, b = r

            xf, yf = local_xy_fractional(
                frame,
                a,
                b,
            )

            alive = residual_representation(
                frame,
                a,
                b,
            )

            print(
                f"    z={z:<2} "
                f"a={a:<6} "
                f"b={b:<6} "
                f"x={str(xf):<8} "
                f"y={str(yf):<8} "
                f"bits=({a & 1},{b & 1}) "
                f"alive={alive}"
            )

            if not alive:
                break

        future = actual_future_steps(
            frame,
            state.p,
            state.q,
            2,
        )

        predicted = None

        r2 = residuals_at_level(
            frame,
            state.p,
            state.q,
            2,
        )

        if r2 is not None:

            a2, b2 = r2

            predicted = predicted_future_steps(
                frame,
                a2,
                b2,
            )

        print(
            f"    future_actual={future}"
        )

        print(
            f"    future_predicted={predicted}"
        )

    print()


# ==============================================================================
# SYMBOLIC SUMMARY
# ==============================================================================

def summary() -> None:

    print("=" * 90)
    print("SYMBOLIC SUMMARY")
    print("=" * 90)

    print(
r"""
At level z let:

    M = 2^(z-1).

FRAME A:

    a = (p-3)/M
    b = (q+3)/M

    x = (b-a)/2
    y = (a+b)/2.

FRAME B:

    a = (p+1)/M
    b = (q-3)/M

    x = (3a-b)/2
    y = (3a+b)/2.

The important distinction is:

    residuals_at_level()
        merely tells us that the factor congruences
        hold at that modulus.

    residual_representation()
        additionally requires x,y to be integers.

Thus:

    a == b (mod 2)

is the current representation condition.

A residual pair such as:

    B: (a,b)=(1,0)

is therefore a legitimate DEAD state, not an exception.

For a live state, the next-level residuals are:

    (a',b')=(a/2,b/2).

The next state is live exactly when:

    (a/2,b/2)

again has equal parity.

The experiment compares this recursive rule with the
directly observed A/B levels and the global invariant
coordinates X,Y.
"""
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 614 START")
    print("=" * 90)
    print()
    print(
        "EXACT 2-ADIC RESIDUAL LADDER "
        "— SAFE INTEGER/DEAD STATE HANDLING"
    )
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

    states = generate_odd_semiprimes(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    print("[3] GLOBAL STATES")

    print(
        f"    states={len(states)}"
    )
    print()

    test_domain(states)
    test_current_representation(states)
    test_one_step(states)
    test_multi_level_ladder(states)
    test_residual_halving(states)
    test_global_crosswalk(states)
    test_depth(states)
    test_minimal_modulus(states)
    test_bit_stream(states)

    examples(states)
    summary()

    print("=" * 90)
    print("EXPERIMENT 614 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()