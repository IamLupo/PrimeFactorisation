#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 615 START
==========================================================================================

CLOSED-FORM FACTOR RESIDUAL DEPTH LAW

Experiment 614 established exactly:

    current representation:
        a == b (mod 2)

    next representation:
        residuals (a/2,b/2)
        must again have equal parity.

The recursive ladder suggests a closed form.

Let:

    alpha = v2(a)
    beta  = v2(b)

For a currently represented residual pair:

    if alpha == beta:
        future_steps = alpha

    if alpha != beta:
        future_steps = min(alpha,beta) - 1

This can be written uniformly as:

    future_steps
        =
    min(
        v2(a),
        v2(b),
        v2(a-b) - 1
    )

with the usual convention:

    v2(0) = +infinity.

Reason:

    if alpha != beta:

        v2(a-b) = min(alpha,beta)

        so:

        future_steps = min(alpha,beta)-1.

    if alpha == beta=t:

        a/2^t and b/2^t are both odd,

        therefore:

        v2(a-b) > t,

        giving:

        future_steps=t.

This experiment tests that law directly.

It also derives the same law from:

    X,Y

and from:

    p,q

for both frames.

==========================================================================================
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


# ==============================================================================
# CONFIG
# ==============================================================================

PRIME_LIMIT = 6250
MAX_LEVEL = 24
MAX_EXAMPLES = 25


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

            sieve[start:limit + 1:p] = (
                b"\x00" * count
            )

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
        f"Unexpected odd-semiprime residue: "
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

        Xn = q - p + 6
        Yn = p + q

    elif frame == "B":

        Xn = 3 * p - q + 6
        Yn = 3 * p + q

    else:
        raise ValueError(frame)

    if Xn & 1 or Yn & 1:

        raise ArithmeticError(
            f"Non-integral global coordinates: "
            f"frame={frame}, p={p}, q={q}"
        )

    return Xn // 2, Yn // 2


# ==============================================================================
# V2
# ==============================================================================

INF = 10**9


def v2(n: int) -> int:

    if n == 0:
        return INF

    n = abs(n)

    return (
        n & -n
    ).bit_length() - 1


def v2_string(n: int) -> str:

    v = v2(n)

    if v >= INF:
        return "inf"

    return str(v)


# ==============================================================================
# RESIDUALS AT LEVEL
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
        raise ValueError(frame)

    if rp % M != 0:
        return None

    if rq % M != 0:
        return None

    return (
        rp // M,
        rq // M,
    )


# ==============================================================================
# LOCAL RESIDUAL COORDINATES
# ==============================================================================

def residual_local_xy(
    frame: str,
    a: int,
    b: int,
) -> tuple[Fraction, Fraction]:

    if frame == "A":

        return (
            Fraction(b - a, 2),
            Fraction(a + b, 2),
        )

    if frame == "B":

        return (
            Fraction(3 * a - b, 2),
            Fraction(3 * a + b, 2),
        )

    raise ValueError(frame)


def residual_is_live(
    frame: str,
    a: int,
    b: int,
) -> bool:

    x, y = residual_local_xy(
        frame,
        a,
        b,
    )

    return (
        x.denominator == 1
        and y.denominator == 1
    )


# ==============================================================================
# DIRECT OBSERVED FUTURE STEPS
# ==============================================================================

def observed_future_steps(
    frame: str,
    p: int,
    q: int,
    z: int,
) -> int:

    steps = 0
    current_z = z

    while current_z < MAX_LEVEL:

        current = residuals_at_level(
            frame,
            p,
            q,
            current_z,
        )

        if current is None:
            break

        a, b = current

        if not residual_is_live(
            frame,
            a,
            b,
        ):
            break

        nxt = residuals_at_level(
            frame,
            p,
            q,
            current_z + 1,
        )

        if nxt is None:
            break

        na, nb = nxt

        if not residual_is_live(
            frame,
            na,
            nb,
        ):
            break

        steps += 1
        current_z += 1

    return steps


# ==============================================================================
# CLOSED FORM RESIDUAL LAW
# ==============================================================================

def residual_future_formula(
    a: int,
    b: int,
) -> int:

    """
    Closed form:

        min(v2(a), v2(b), v2(a-b)-1)

    provided (a,b) is currently representable.

    The caller must ensure:

        a == b mod 2.
    """

    if (a - b) & 1:

        return -1

    va = v2(a)
    vb = v2(b)
    vd = v2(a - b)

    terms = [
        va,
        vb,
    ]

    if vd >= INF:
        # a == b.
        #
        # In this case v2(a-b)=inf, so it does not
        # constrain the minimum.
        pass

    else:
        terms.append(vd - 1)

    return min(terms)


# ==============================================================================
# GLOBAL DEPTH
# ==============================================================================

def global_depth(
    X: int,
    Y: int,
) -> int:

    return (
        1
        + min(
            v2(X),
            v2(Y),
        )
    )


def observed_depth(
    frame: str,
    p: int,
    q: int,
) -> int:

    return (
        2
        + observed_future_steps(
            frame,
            p,
            q,
            2,
        )
    )


# ==============================================================================
# FACTOR-LEVEL CLOSED FORM DIRECTLY FROM p,q
# ==============================================================================

def base_residuals(
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

    raise ValueError(frame)


def factor_depth_formula(
    frame: str,
    p: int,
    q: int,
) -> int:

    a, b = base_residuals(
        frame,
        p,
        q,
    )

    future = residual_future_formula(
        a,
        b,
    )

    if future < 0:
        raise ArithmeticError(
            f"Base residual pair not representable: "
            f"frame={frame}, a={a}, b={b}"
        )

    return 2 + future


# ==============================================================================
# TEST 1
# ==============================================================================

def test_base_formula(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 1: BASE RESIDUAL CLOSED FORM")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        a, b = base_residuals(
            frame,
            state.p,
            state.q,
        )

        if not residual_is_live(
            frame,
            a,
            b,
        ):
            failures += 1

            if failures <= MAX_EXAMPLES:
                print(
                    f"    unexpected dead base state "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"a={a} b={b}"
                )

            continue

        actual = observed_future_steps(
            frame,
            state.p,
            state.q,
            2,
        )

        predicted = residual_future_formula(
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
                    f"frame={frame} "
                    f"a={a} b={b} "
                    f"actual={actual} "
                    f"predicted={predicted} "
                    f"v2a={v2_string(a)} "
                    f"v2b={v2_string(b)} "
                    f"v2diff={v2_string(a-b)}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 2
# ==============================================================================

def test_current_level_formula(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 2: CURRENT-LEVEL CLOSED FORM")
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

            residuals = residuals_at_level(
                frame,
                state.p,
                state.q,
                z,
            )

            if residuals is None:
                continue

            a, b = residuals

            if not residual_is_live(
                frame,
                a,
                b,
            ):
                continue

            actual = observed_future_steps(
                frame,
                state.p,
                state.q,
                z,
            )

            predicted = residual_future_formula(
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
                        f"a={a} b={b} "
                        f"actual={actual} "
                        f"predicted={predicted} "
                        f"v2a={v2_string(a)} "
                        f"v2b={v2_string(b)} "
                        f"v2diff={v2_string(a-b)}"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 3
# ==============================================================================

def test_depth_formula(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 3: FACTOR DEPTH == GLOBAL DEPTH")
    print("=" * 90)

    checked = 0
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

        actual = observed_depth(
            frame,
            state.p,
            state.q,
        )

        global_pred = global_depth(
            X,
            Y,
        )

        factor_pred = factor_depth_formula(
            frame,
            state.p,
            state.q,
        )

        distribution[actual] += 1

        checked += 1

        if actual != global_pred:

            failures += 1

            if failures <= MAX_EXAMPLES:
                print(
                    f"    global mismatch "
                    f"n={state.n} "
                    f"actual={actual} "
                    f"global={global_pred}"
                )

        if actual != factor_pred:

            failures += 1

            if failures <= MAX_EXAMPLES:
                print(
                    f"    factor mismatch "
                    f"n={state.n} "
                    f"actual={actual} "
                    f"factor={factor_pred}"
                )

    print()
    print("DEPTH DISTRIBUTION")

    for depth in sorted(
        distribution
    ):
        print(
            f"    depth={depth}: "
            f"{distribution[depth]}"
        )

    print()

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 4
# ==============================================================================

def test_valuation_cases(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 4: EQUAL VS UNEQUAL v2 CASES")
    print("=" * 90)

    checked = 0
    failures = 0

    case_counts = Counter()

    for state in states:

        frame = frame_from_n(
            state.n
        )

        a, b = base_residuals(
            frame,
            state.p,
            state.q,
        )

        va = v2(a)
        vb = v2(b)

        if va == vb:
            case = "equal"
            expected = va
        else:
            case = "unequal"

            expected = (
                min(va, vb) - 1
            )

        actual = residual_future_formula(
            a,
            b,
        )

        checked += 1
        case_counts[case] += 1

        if actual != expected:

            failures += 1

            if failures <= MAX_EXAMPLES:
                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"a={a} b={b} "
                    f"va={v2_string(a)} "
                    f"vb={v2_string(b)} "
                    f"actual={actual} "
                    f"expected={expected}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print(
        f"    equal valuation states="
        f"{case_counts['equal']}"
    )

    print(
        f"    unequal valuation states="
        f"{case_counts['unequal']}"
    )

    print()


# ==============================================================================
# TEST 5
# ==============================================================================

def test_gcd_version(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 5: gcd / 2-ADIC VERSION")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        a, b = base_residuals(
            frame,
            state.p,
            state.q,
        )

        if not residual_is_live(
            frame,
            a,
            b,
        ):
            continue

        future = residual_future_formula(
            a,
            b,
        )

        g = abs(
            __import__("math").gcd(
                a,
                b,
                a - b,
            )
        )

        vg = v2(g)

        # For a currently live pair:
        #
        # if both valuations equal t:
        #     v2(g)=t
        #     future=t
        #
        # if unequal:
        #     v2(g)=min(v2(a),v2(b))
        #     future=v2(g)-1.
        #
        # Hence future differs from v2(g) by exactly one in
        # the unequal case.

        va = v2(a)
        vb = v2(b)

        expected = (
            vg
            if va == vb
            else vg - 1
        )

        checked += 1

        if future != expected:

            failures += 1

            if failures <= MAX_EXAMPLES:
                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"a={a} b={b} "
                    f"g={g} "
                    f"v2g={v2_string(g)} "
                    f"va={v2_string(a)} "
                    f"vb={v2_string(b)} "
                    f"future={future} "
                    f"expected={expected}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 6
# ==============================================================================

def test_recursive_shift(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 6: RESIDUAL SHIFT / RECURSIVE VALUATION")
    print("=" * 90)

    checked = 0
    failures = 0

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

            if not residual_is_live(
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
                        f"r=({a},{b}) "
                        f"next=({ap},{bp}) "
                        f"expected=({a//2},{b//2})"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 7
# ==============================================================================

def test_depth_from_pq(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 7: DIRECT p,q DEPTH FORMULA")
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        factor_depth = factor_depth_formula(
            frame,
            state.p,
            state.q,
        )

        actual = observed_depth(
            frame,
            state.p,
            state.q,
        )

        checked += 1

        if factor_depth != actual:

            failures += 1

            if failures <= MAX_EXAMPLES:
                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={frame} "
                    f"p={state.p} "
                    f"q={state.q} "
                    f"factor={factor_depth} "
                    f"actual={actual}"
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
        s.n: s
        for s in states
        if s.n in wanted
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

        a, b = base_residuals(
            frame,
            state.p,
            state.q,
        )

        future = residual_future_formula(
            a,
            b,
        )

        depth = 2 + future

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
            f"    v2X={v2_string(X)} "
            f"v2Y={v2_string(Y)}"
        )

        print(
            f"    base residuals="
            f"({a},{b})"
        )

        print(
            f"    v2(a)={v2_string(a)} "
            f"v2(b)={v2_string(b)} "
            f"v2(a-b)={v2_string(a-b)}"
        )

        print(
            f"    future_steps={future}"
        )

        print(
            f"    predicted_depth={depth}"
        )

        print()

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

            az, bz = r

            x, y = residual_local_xy(
                frame,
                az,
                bz,
            )

            alive = residual_is_live(
                frame,
                az,
                bz,
            )

            next_alive = False

            if alive and z < MAX_LEVEL:

                rn = residuals_at_level(
                    frame,
                    state.p,
                    state.q,
                    z + 1,
                )

                if rn is not None:

                    na, nb = rn

                    next_alive = residual_is_live(
                        frame,
                        na,
                        nb,
                    )

            print(
                f"    z={z:<2} "
                f"a={az:<7} "
                f"b={bz:<7} "
                f"x={str(x):<8} "
                f"y={str(y):<8} "
                f"bits=({az & 1},{bz & 1}) "
                f"alive={alive} "
                f"next={next_alive}"
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
At level z define:

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

CURRENT REPRESENTATION

    A or B is integer-representable exactly when:

        a == b (mod 2).

NEXT LEVEL

    residuals become:

        a' = a/2
        b' = b/2

    so another level survives exactly when:

        a/2 == b/2 (mod 2).

CLOSED FORM

Let:

    alpha = v2(a)
    beta  = v2(b)
    gamma = v2(a-b).

For a currently represented state:

    future_steps
        =
    min(alpha, beta, gamma-1).

Equivalent piecewise form:

    if alpha == beta:

        future_steps = alpha

    if alpha != beta:

        future_steps = min(alpha,beta)-1.

This explains:

    n=21:
        (a,b)=(2,2)
        alpha=beta=1
        future=1

    n=69:
        (a,b)=(2,10)
        alpha=beta=1
        future=1

    n=57:
        (a,b)=(2,8)
        alpha=1, beta=3
        future=0

    n=15:
        (a,b)=(0,4)
        alpha=inf, beta=2
        future=1.

GLOBAL DEPTH

At the base level z=2:

    depth = 2 + future_steps.

This must equal:

    depth = 1 + min(v2(X),v2(Y)).

Thus the experiment tests the full equality:

    2
    + min(
        v2(p-residue),
        v2(q-residue),
        v2(difference-residual)-1
      )

        =

    1 + min(v2(X),v2(Y)).

The important result would be a direct closed form
for the entire hierarchy using only valuations of
simple linear forms in the prime factors.
"""
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 615 START")
    print("=" * 90)
    print()
    print(
        "CLOSED-FORM FACTOR RESIDUAL DEPTH LAW"
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

    test_base_formula(states)
    test_current_level_formula(states)
    test_depth_formula(states)
    test_valuation_cases(states)
    test_gcd_version(states)
    test_recursive_shift(states)
    test_depth_from_pq(states)

    examples(states)
    symbolic_summary()

    print("=" * 90)
    print("EXPERIMENT 615 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
