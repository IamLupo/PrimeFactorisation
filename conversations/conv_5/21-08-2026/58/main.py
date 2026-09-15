#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 619
# ==============================================================================
#
# COLLAPSE THE EQUALITY FLAG
#
# Experiment 617:
#
#     depth = min(v2(A),v2(B)) + [v2(A)=v2(B)]
#
# Experiment 618:
#
#     [v2(A)=v2(B)]
#
# is exactly equivalent to:
#
#     v2(A+B) > min(v2(A),v2(B))
#
# and also:
#
#     v2(A-B) > min(v2(A),v2(B)).
#
# This experiment searches for an even smaller formulation.
#
# FRAME A:
#
#     A = p-3
#     B = q+3
#
#     A+B = p+q
#     A-B = p-q-6
#
# FRAME B:
#
#     A = p+1
#     B = q-3
#
#     A+B = p+q-2
#     A-B = p-q+4
#
# Goals:
#
#   1. Determine whether E can be recovered from one valuation
#      of a simple polynomial.
#
#   2. Determine whether E can be recovered from gcd(A,B).
#
#   3. Determine whether E can be recovered from p+q or p-q
#      at a modulus determined only by min(v2(A),v2(B)).
#
#   4. Determine whether E is visible directly in n modulo a
#      small state-dependent modulus.
#
#   5. Find the smallest exact Boolean representation.
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250
MAX_BITS = 14
EXAMPLE_NS = [
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
    141,
    183,
    213,
]


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


INF = 10**9


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):

        if sieve[p]:

            start = p * p
            count = ((limit - start) // p) + 1

            sieve[
                start:limit + 1:p
            ] = b"\x00" * count

    return [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_states(
    primes: list[int],
) -> list[State]:

    states = []

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
        f"Unexpected odd semiprime residue: "
        f"n={n}, n mod 4={r}"
    )


# ==============================================================================
# RAW RESIDUALS
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

    raise ValueError(frame)


# ==============================================================================
# V2
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return INF

    x = abs(x)

    return (
        x & -x
    ).bit_length() - 1


def v2g(a: int, b: int) -> int:

    return v2(
        gcd(
            abs(a),
            abs(b),
        )
    )


# ==============================================================================
# ODD PART
# ==============================================================================

def odd_part(x: int) -> int:

    if x == 0:
        return 0

    return abs(x) >> v2(x)


# ==============================================================================
# EQUALITY FLAG
# ==============================================================================

def equality_flag(
    A: int,
    B: int,
) -> bool:

    return v2(A) == v2(B)


# ==============================================================================
# DEPTH
# ==============================================================================

def depth_formula(
    A: int,
    B: int,
) -> int:

    alpha = v2(A)
    beta = v2(B)

    return (
        min(alpha, beta)
        + int(alpha == beta)
    )


# ==============================================================================
# FORMULAS
# ==============================================================================

def derived_forms(
    frame: str,
    p: int,
    q: int,
) -> dict[str, int]:

    A, B = raw_residuals(
        frame,
        p,
        q,
    )

    return {
        "A": A,
        "B": B,
        "A+B": A + B,
        "A-B": A - B,
        "B-A": B - A,
        "p+q": p + q,
        "p-q": p - q,
        "q-p": q - p,
        "n": p * q,
    }


# ==============================================================================
# TEST 1
# ==============================================================================

def test_exact_sum_difference(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 1: SUM / DIFFERENCE COLLAPSE")
    print("=" * 90)

    candidates = {
        "v2(A+B)>min": lambda A, B:
            v2(A + B) > min(v2(A), v2(B)),

        "v2(A-B)>min": lambda A, B:
            v2(A - B) > min(v2(A), v2(B)),

        "v2(A+B)==min+1": lambda A, B:
            v2(A + B) == min(v2(A), v2(B)) + 1,

        "v2(A-B)==min+1": lambda A, B:
            v2(A - B) == min(v2(A), v2(B)) + 1,
    }

    for name, predicate in candidates.items():

        failures = 0

        for state in states:

            frame = frame_from_n(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            actual = equality_flag(
                A,
                B,
            )

            predicted = predicate(
                A,
                B,
            )

            if actual != predicted:
                failures += 1

        print(
            f"    {name:<28} "
            f"failures={failures}"
        )

    print()


# ==============================================================================
# TEST 2
# ==============================================================================

def test_linear_factor_forms(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 2: SINGLE LINEAR-FORM COLLAPSE")
    print("=" * 90)

    # A simple single valuation predicate can only use one
    # integer-valued form plus the minimum valuation.
    #
    # Search:
    #
    #     v2(form) > min
    #
    # for several equivalent linear forms.

    for name in [
        "A+B",
        "A-B",
        "B-A",
        "p+q",
        "p-q",
        "q-p",
    ]:

        failures = 0

        for state in states:

            frame = frame_from_n(
                state.n
            )

            forms = derived_forms(
                frame,
                state.p,
                state.q,
            )

            A = forms["A"]
            B = forms["B"]

            actual = equality_flag(
                A,
                B,
            )

            predicted = (
                v2(
                    forms[name]
                )
                >
                min(
                    v2(A),
                    v2(B),
                )
            )

            if actual != predicted:
                failures += 1

        print(
            f"    {name:<8} "
            f"failures={failures}"
        )

    print()


# ==============================================================================
# TEST 3
# ==============================================================================

def test_gcd_equivalence(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 3: GCD COLLAPSE")
    print("=" * 90)

    candidates = {
        "gcd(A,B) has exact min":
            lambda A, B:
                v2g(A, B)
                ==
                min(v2(A), v2(B)),

        "A/g and B/g odd":
            lambda A, B:
                (
                    A == 0
                    or B == 0
                    or (
                        (A // gcd(abs(A), abs(B))) & 1
                        and
                        (B // gcd(abs(A), abs(B))) & 1
                    )
                ),

        "A/g and B/g both odd exactly":
            lambda A, B:
                (
                    A != 0
                    and B != 0
                    and
                    ((A // gcd(abs(A), abs(B))) & 1) == 1
                    and
                    ((B // gcd(abs(A), abs(B))) & 1) == 1
                ),
    }

    for name, predicate in candidates.items():

        failures = 0

        for state in states:

            frame = frame_from_n(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            actual = equality_flag(
                A,
                B,
            )

            predicted = predicate(
                A,
                B,
            )

            if actual != predicted:
                failures += 1

        print(
            f"    {name:<35} "
            f"failures={failures}"
        )

    print()


# ==============================================================================
# TEST 4
# ==============================================================================

def test_state_dependent_congruence(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 4: STATE-DEPENDENT CONGRUENCE")
    print("=" * 90)

    # Let:
    #
    #     m = min(v2(A),v2(B)).
    #
    # Equality means:
    #
    #     A == 0 mod 2^(m+1)
    #     B == 0 mod 2^(m+1)
    #
    # but unequal valuations fail one of these.
    #
    # Check several equivalent formulations.

    failures_1 = 0
    failures_2 = 0
    failures_3 = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        m = min(
            v2(A),
            v2(B),
        )

        actual = equality_flag(
            A,
            B,
        )

        if m >= INF:
            predicted_1 = True
        else:
            modulus = 1 << (m + 1)

            predicted_1 = (
                A % modulus == 0
                and
                B % modulus == 0
            )

            predicted_2 = (
                (A + B) % modulus == 0
                and
                (A - B) % modulus == 0
            )

            predicted_3 = (
                (A + B) % modulus == 0
            )

            if actual != predicted_2:
                failures_2 += 1

            if actual != predicted_3:
                failures_3 += 1

        if actual != predicted_1:
            failures_1 += 1

    print(
        "    A,B divisible by 2^(m+1)"
    )
    print(
        f"        failures={failures_1}"
    )

    print(
        "    A+B and A-B divisible by 2^(m+1)"
    )
    print(
        f"        failures={failures_2}"
    )

    print(
        "    A+B divisible by 2^(m+1)"
    )
    print(
        f"        failures={failures_3}"
    )

    print()


# ==============================================================================
# TEST 5
# ==============================================================================

def test_n_congruence(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 5: CAN EQUALITY BE READ DIRECTLY FROM n?")
    print("=" * 90)

    # At each state we know:
    #
    #     m = min(v2(A),v2(B)).
    #
    # Search whether E is determined by n modulo:
    #
    #     2^(m+1)
    #     2^(m+2)
    #     2^(m+3)
    #
    # This is state-dependent, but no longer explicitly uses
    # the two residuals beyond m.

    for extra in range(1, 7):

        table = defaultdict(set)

        for state in states:

            frame = frame_from_n(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            m = min(
                v2(A),
                v2(B),
            )

            if m >= INF:
                continue

            modulus = 1 << (
                m + extra
            )

            key = (
                frame,
                m,
                state.n % modulus,
            )

            table[key].add(
                equality_flag(A, B)
            )

        ambiguous = sum(
            len(values) > 1
            for values in table.values()
        )

        print(
            f"    extra={extra} "
            f"modulus=2^(m+{extra}) "
            f"states={len(table)} "
            f"ambiguous={ambiguous}"
        )

    print()


# ==============================================================================
# TEST 6
# ==============================================================================

def test_pq_congruence(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 6: FACTOR CONGRUENCE COLLAPSE")
    print("=" * 90)

    # Search whether equality is exactly determined by:
    #
    #     p+q mod 2^(m+1)
    #
    # or:
    #
    #     p-q + frame_constant mod 2^(m+1).
    #
    # The frame is retained because the two coordinate systems
    # have different offsets.

    for extra in range(1, 6):

        failures_sum = 0
        failures_diff = 0

        for state in states:

            frame = frame_from_n(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            m = min(
                v2(A),
                v2(B),
            )

            if m >= INF:
                continue

            modulus = 1 << (
                m + extra
            )

            actual = equality_flag(
                A,
                B,
            )

            if frame == "A":

                plus_form = (
                    state.p + state.q
                )

                minus_form = (
                    state.p
                    - state.q
                    - 6
                )

            else:

                plus_form = (
                    state.p
                    + state.q
                    - 2
                )

                minus_form = (
                    state.p
                    - state.q
                    + 4
                )

            # Equality implies the relevant residual forms
            # vanish one bit further than m.
            #
            # For extra=1 this is the sharp condition.
            #
            pred_sum = (
                plus_form % modulus == 0
            )

            pred_diff = (
                minus_form % modulus == 0
            )

            if pred_sum != actual:
                failures_sum += 1

            if pred_diff != actual:
                failures_diff += 1

        print(
            f"    extra={extra} "
            f"sum_failures={failures_sum} "
            f"diff_failures={failures_diff}"
        )

    print()


# ==============================================================================
# TEST 7
# ==============================================================================

def test_xy_equivalence(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 7: EQUALITY FLAG THROUGH X,Y")
    print("=" * 90)

    # Derive X,Y directly and search whether equality corresponds
    # to a simple relation between v2(X), v2(Y) and odd parts.
    #
    # Since:
    #
    #     depth = 1 + min(v2(X),v2(Y))
    #
    # while:
    #
    #     depth = min(alpha,beta) + E,
    #
    # equality must occur exactly when:
    #
    #     min(v2(X),v2(Y))
    #       =
    #     min(alpha,beta)
    #
    # This test explicitly verifies that relation.

    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)

        X, Y = (
            (
                state.q - state.p + 6
            ) // 2,
            (
                state.p + state.q
            ) // 2,
        ) if frame == "A" else (
            (
                3 * state.p
                - state.q
                + 6
            ) // 2,
            (
                3 * state.p
                + state.q
            ) // 2,
        )

        E = (
            alpha == beta
        )

        predicted_depth = (
            min(alpha, beta)
            + int(E)
        )

        global_depth = (
            1
            + min(
                v2(X),
                v2(Y),
            )
        )

        if predicted_depth != global_depth:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 8
# ==============================================================================

def test_candidate_boolean_normal_forms(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 8: BOOLEAN NORMAL-FORM SEARCH")
    print("=" * 90)

    # Build a small feature basis:
    #
    #   s1 = [v2(A+B) > m]
    #   s2 = [v2(A-B) > m]
    #   s3 = [(A/g)%2 == 1]
    #   s4 = [(B/g)%2 == 1]
    #
    # Then test all Boolean functions of these four features
    # against the equality flag. This is intentionally exhaustive
    # over the 2^16 truth tables.
    #
    # This tests whether some surprisingly simple low-level
    # Boolean expression exists.
    #
    # In practice s1 and s2 should already equal E, but the
    # enumeration also verifies redundancy.

    rows = []

    for state in states:

        frame = frame_from_n(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        m = min(
            v2(A),
            v2(B),
        )

        s1 = (
            v2(A + B) > m
        )

        s2 = (
            v2(A - B) > m
        )

        g = gcd(
            abs(A),
            abs(B),
        )

        if g == 0:

            s3 = False
            s4 = False

        else:

            s3 = (
                A % (2 * g) == 0
            )

            s4 = (
                B % (2 * g) == 0
            )

        target = equality_flag(
            A,
            B,
        )

        rows.append(
            (
                (
                    int(s1),
                    int(s2),
                    int(s3),
                    int(s4),
                ),
                target,
            )
        )

    pattern = {}

    ambiguous = set()

    for features, target in rows:

        if features not in pattern:

            pattern[features] = target

        elif pattern[features] != target:

            ambiguous.add(features)

    print(
        f"feature states={len(pattern)} "
        f"ambiguous={len(ambiguous)}"
    )

    for key in sorted(pattern):

        print(
            f"    {key} -> {pattern[key]}"
        )

    print()


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states: list[State],
) -> None:

    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    lookup = {
        state.n: state
        for state in states
    }

    for n in EXAMPLE_NS:

        state = lookup.get(n)

        if state is None:
            continue

        frame = frame_from_n(
            n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)
        m = min(alpha, beta)

        if frame == "A":

            plus_form = (
                state.p + state.q
            )

            diff_form = (
                state.p
                - state.q
                - 6
            )

        else:

            plus_form = (
                state.p
                + state.q
                - 2
            )

            diff_form = (
                state.p
                - state.q
                + 4
            )

        equal = (
            alpha == beta
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
            f"    A={A} "
            f"B={B}"
        )

        print(
            f"    v2(A)={alpha if alpha < INF else 'inf'} "
            f"v2(B)={beta if beta < INF else 'inf'}"
        )

        print(
            f"    min={m if m < INF else 'inf'} "
            f"equal={equal}"
        )

        print(
            f"    A+B={plus_form}"
        )

        print(
            f"    A-B={diff_form}"
        )

        print(
            f"    v2(A+B)="
            f"{v2(plus_form) if v2(plus_form) < INF else 'inf'}"
        )

        print(
            f"    v2(A-B)="
            f"{v2(diff_form) if v2(diff_form) < INF else 'inf'}"
        )

        if m < INF:

            modulus = 1 << (
                m + 1
            )

            print(
                f"    2^(m+1)={modulus}"
            )

            print(
                f"    A mod 2^(m+1)="
                f"{A % modulus}"
            )

            print(
                f"    B mod 2^(m+1)="
                f"{B % modulus}"
            )

            print(
                f"    (A+B) mod="
                f"{plus_form % modulus}"
            )

            print(
                f"    (A-B) mod="
                f"{diff_form % modulus}"
            )

        print(
            f"    depth="
            f"{depth_formula(A,B)}"
        )

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
The current exact law is:

    depth =
        min(alpha,beta)
        + [alpha=beta]

where:

    alpha=v2(A)
    beta =v2(B).

For unequal valuations:

    alpha != beta

    =>

    v2(A+B)
      =
    v2(A-B)
      =
    min(alpha,beta).

For equal valuations:

    alpha=beta=t

    write:

        A=2^t a
        B=2^t b

    with a,b odd.

    Therefore:

        A+B=2^t(a+b)
        A-B=2^t(a-b)

    and both a+b and a-b are even.

Hence:

    [alpha=beta]

        <=>

    v2(A+B)>min(alpha,beta)

        <=>

    v2(A-B)>min(alpha,beta).

For the two frames:

    FRAME A:

        A=p-3
        B=q+3

        A+B=p+q
        A-B=p-q-6

    FRAME B:

        A=p+1
        B=q-3

        A+B=p+q-2
        A-B=p-q+4.

So the target collapse is:

    depth =
        min(v2(A),v2(B))
        +
        [
            v2(L)>min(v2(A),v2(B))
        ]

for some single simple linear form L.

If Experiment 619 finds:

    L = p+q
    or
    L = p+q-2
    or
    L = p-q-6
    or
    L = p-q+4,

with the appropriate frame-dependent substitution,
then the equality branch is completely reducible to
a single extra 2-adic test.

An even stronger result would be an n-only formulation:

    E =
    predicate(
        n mod 2^(m+c),
        m,
        frame
    )

which would eliminate explicit p and q from the depth logic.
"""
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 619 START")
    print("=" * 90)
    print()
    print("COLLAPSE THE EQUALITY FLAG")
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
        f"    states={len(states)}"
    )
    print()

    test_exact_sum_difference(
        states
    )

    test_linear_factor_forms(
        states
    )

    test_gcd_equivalence(
        states
    )

    test_state_dependent_congruence(
        states
    )

    test_n_congruence(
        states
    )

    test_pq_congruence(
        states
    )

    test_xy_equivalence(
        states
    )

    test_candidate_boolean_normal_forms(
        states
    )

    print_examples(
        states
    )

    symbolic_summary()

    print("=" * 90)
    print("EXPERIMENT 619 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
