#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 618
# ==============================================================================
#
# EQUALITY-BRANCH COLLAPSE
#
# Experiment 617 established the exact depth:
#
#     depth =
#         min(v2(A), v2(B))
#         + [v2(A) == v2(B)].
#
# where:
#
# FRAME A:
#     A = p - 3
#     B = q + 3
#
# FRAME B:
#     A = p + 1
#     B = q - 3
#
# The only remaining nonlinear-looking component is:
#
#     E = [v2(A) == v2(B)].
#
# This experiment studies that equality flag.
#
# For equal valuations:
#
#     A = 2^t * a
#     B = 2^t * b
#
# with a,b odd.
#
# We test whether E can be characterized by:
#
#     a+b mod 4
#     a-b mod 4
#     a*b mod 4
#     a+b mod 8
#     a-b mod 8
#     p,q modulo small powers of 2
#     p/q relationships
#     gcd(A,B)
#     gcd(A+B,A-B)
#     X,Y valuations
#     X/Y valuation equality
#
# We also test a fundamental identity:
#
#     A+B
#     A-B
#
# because:
#
#     v2(A)=v2(B)
#
# may have a simple interpretation through the valuations of
#     A+B and A-B.
#
# For equal alpha=beta=t:
#
#     v2(A+B) >= t+1
#     v2(A-B) >= t+1.
#
# For unequal alpha != beta:
#
#     min(v2(A+B), v2(A-B)) = min(alpha,beta).
#
# The experiment searches for the smallest exact predicate.
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250
MAX_MODULUS_POWER = 12
EXAMPLE_LIMIT = 30


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

    if limit < 3:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):

        if sieve[p]:

            start = p * p
            count = ((limit - start) // p) + 1

            sieve[
                start : limit + 1 : p
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
        f"Expected odd semiprime, "
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

        xn = q - p + 6
        yn = p + q

    elif frame == "B":

        xn = 3 * p - q + 6
        yn = 3 * p + q

    else:
        raise ValueError(frame)

    if xn % 2 or yn % 2:

        raise ArithmeticError(
            f"Non-integral X,Y: "
            f"frame={frame} "
            f"p={p} "
            f"q={q}"
        )

    return xn // 2, yn // 2


# ==============================================================================
# RAW RESIDUALS
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

def v2(x: int) -> int:

    if x == 0:
        return INF

    x = abs(x)

    return (
        x & -x
    ).bit_length() - 1


def vstr(x: int) -> str:

    value = v2(x)

    if value >= INF:
        return "inf"

    return str(value)


# ==============================================================================
# ODD PART
# ==============================================================================

def odd_part(x: int) -> int:

    if x == 0:
        return 0

    value = abs(x)
    value >>= v2(value)

    return value


# ==============================================================================
# MOD NORMALIZATION
# ==============================================================================

def mod_pair(
    a: int,
    b: int,
    modulus: int,
) -> tuple[int, int]:

    return (
        a % modulus,
        b % modulus,
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_equality_flag_from_basic_data(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 1: EQUALITY FLAG BASELINE")
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

        equality = (
            alpha == beta
        )

        # Baseline is the definition.
        predicted = (
            v2(A) == v2(B)
        )

        checked += 1

        if equality != predicted:
            failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 2
# ==============================================================================

def test_sum_difference_structure(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 2: SUM / DIFFERENCE 2-ADIC STRUCTURE")
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

        equal = (
            alpha == beta
        )

        if equal:

            gamma_sum = v2(A + B)
            gamma_diff = v2(A - B)

            # For A=2^t*a, B=2^t*b with a,b odd:
            #
            # a+b and a-b are both even.
            #
            # Therefore both valuations must exceed t.
            #
            t = alpha

            if gamma_sum <= t:
                failures += 1

            if gamma_diff <= t:
                failures += 1

        else:

            m = min(alpha, beta)

            gamma_sum = v2(A + B)
            gamma_diff = v2(A - B)

            # With unequal valuations, both A+B and A-B
            # retain the smaller 2-adic valuation.
            #
            if gamma_sum != m:
                failures += 1

            if gamma_diff != m:
                failures += 1

        checked += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 3
# ==============================================================================

def test_odd_part_congruences(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 3: ODD-PART CONGRUENCE SEARCH")
    print("=" * 90)

    # candidate predicates:
    #
    #   a+b mod m
    #   a-b mod m
    #   a*b mod m
    #
    # For equal valuation states only.

    moduli = [
        4,
        8,
        16,
        32,
        64,
        128,
    ]

    predicates = {}

    for modulus in moduli:

        predicates[
            f"a+b mod {modulus}"
        ] = defaultdict(set)

        predicates[
            f"a-b mod {modulus}"
        ] = defaultdict(set)

        predicates[
            f"a*b mod {modulus}"
        ] = defaultdict(set)

    equality_rows = []

    for state in states:

        frame = frame_from_n(state.n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)

        if alpha != beta:
            continue

        a = odd_part(A)
        b = odd_part(B)

        row = {
            "frame": frame,
            "a": a,
            "b": b,
        }

        equality_rows.append(row)

        for modulus in moduli:

            predicates[
                f"a+b mod {modulus}"
            ][
                (a + b) % modulus
            ].add(True)

            predicates[
                f"a-b mod {modulus}"
            ][
                (a - b) % modulus
            ].add(True)

            predicates[
                f"a*b mod {modulus}"
            ][
                (a * b) % modulus
            ].add(True)

    print(
        f"equal-valuation states="
        f"{len(equality_rows)}"
    )

    print()

    for name, mapping in predicates.items():

        residues = sorted(
            mapping.keys()
        )

        print(
            f"    {name}: "
            f"{residues}"
        )

    print()


# ==============================================================================
# TEST 4
# ==============================================================================

def test_small_modulus_state(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 4: CAN EQUALITY BE DETERMINED FROM LOW BITS?")
    print("=" * 90)

    # Compare:
    #
    #   normalized odd parts mod 2^k
    #
    # with equality of valuations.
    #
    # We search for the first k where:
    #
    #   key(Aodd mod 2^k, Bodd mod 2^k)
    #
    # deterministically identifies equal/unequal.
    #
    # Note that the valuation level itself is deliberately excluded;
    # this tests whether odd residual structure contains extra information.

    for k in range(
        1,
        MAX_MODULUS_POWER + 1,
    ):

        modulus = 1 << k

        table = defaultdict(set)

        for state in states:

            frame = frame_from_n(state.n)

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            alpha = v2(A)
            beta = v2(B)

            if alpha >= INF:
                # Treat zero specially.
                odd_a = 0
            else:
                odd_a = odd_part(A)

            if beta >= INF:
                odd_b = 0
            else:
                odd_b = odd_part(B)

            key = (
                odd_a % modulus,
                odd_b % modulus,
            )

            table[key].add(
                alpha == beta
            )

        ambiguous = sum(
            1
            for values in table.values()
            if len(values) > 1
        )

        print(
            f"bits={k:<2} "
            f"modulus={modulus:<5} "
            f"states={len(table):<6} "
            f"ambiguous={ambiguous}"
        )

    print()


# ==============================================================================
# TEST 5
# ==============================================================================

def test_equality_vs_XY(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 5: EQUALITY FLAG VS GLOBAL X,Y")
    print("=" * 90)

    table = defaultdict(set)

    for state in states:

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

        equality = (
            v2(A) == v2(B)
        )

        key = (
            v2(X),
            v2(Y),
            (abs(X) & 7),
            (abs(Y) & 7),
        )

        table[key].add(
            equality
        )

    ambiguous = sum(
        1
        for values in table.values()
        if len(values) > 1
    )

    print(
        f"states={len(table)} "
        f"ambiguous={ambiguous}"
    )

    print()


# ==============================================================================
# TEST 6
# ==============================================================================

def test_equality_vs_prime_residues(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 6: EQUALITY FLAG FROM p,q LOW BITS")
    print("=" * 90)

    for k in range(
        2,
        MAX_MODULUS_POWER + 1,
    ):

        modulus = 1 << k

        table = defaultdict(set)

        for state in states:

            frame = frame_from_n(state.n)

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            key = (
                state.p % modulus,
                state.q % modulus,
                frame,
            )

            table[key].add(
                v2(A) == v2(B)
            )

        ambiguous = sum(
            1
            for values in table.values()
            if len(values) > 1
        )

        print(
            f"bits={k:<2} "
            f"mod={modulus:<5} "
            f"states={len(table):<7} "
            f"ambiguous={ambiguous}"
        )

    print()


# ==============================================================================
# TEST 7
# ==============================================================================

def test_gcd_structure(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 7: GCD STRUCTURE OF EQUALITY BRANCH")
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

        g = gcd(
            abs(A),
            abs(B),
        )

        gval = v2(g)

        checked += 1

        if gval != min(alpha, beta):
            failures += 1

        if alpha == beta:

            # Equal valuation means after removing the gcd's
            # 2-adic part, both normalized quantities are odd.
            #
            if alpha < INF and beta < INF:

                if (
                    (A >> alpha) & 1
                ) != 1:

                    failures += 1

                if (
                    (B >> beta) & 1
                ) != 1:

                    failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 8
# ==============================================================================

def test_direct_equality_formula(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 8: SEARCH FOR SIMPLE EQUALITY PREDICATES")
    print("=" * 90)

    # Candidate exact predicates.
    #
    # Each candidate maps a state to True/False and is compared
    # to alpha == beta.
    #
    # The point is not merely to validate the definition but to
    # identify alternate forms that can remove the explicit
    # valuation comparison.

    candidates = {
        "v2(A)==v2(B)": lambda A, B:
            v2(A) == v2(B),

        "v2(A+B)>min": lambda A, B:
            v2(A + B) > min(v2(A), v2(B)),

        "v2(A-B)>min": lambda A, B:
            v2(A - B) > min(v2(A), v2(B)),

        "both_sum_and_diff_gt_min":
            lambda A, B:
                (
                    v2(A + B) > min(v2(A), v2(B))
                    and
                    v2(A - B) > min(v2(A), v2(B))
                ),

        "v2(gcd)=v2(A)=v2(B)":
            lambda A, B:
                (
                    v2(gcd(abs(A), abs(B)))
                    == v2(A)
                    == v2(B)
                ),
    }

    for name, predicate in candidates.items():

        failures = 0

        for state in states:

            frame = frame_from_n(state.n)

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            actual = (
                v2(A) == v2(B)
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
        141,
        183,
        213,
    ]

    lookup = {
        s.n: s
        for s in states
    }

    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    for n in requested:

        state = lookup.get(n)

        if state is None:
            continue

        frame = frame_from_n(n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)

        equal = (
            alpha == beta
        )

        gcd_val = v2(
            gcd(
                abs(A),
                abs(B),
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
            f"    A={A} "
            f"B={B}"
        )

        print(
            f"    v2(A)={vstr(A)} "
            f"v2(B)={vstr(B)} "
            f"equal={equal}"
        )

        if A != 0:
            print(
                f"    odd(A)={odd_part(A)}"
            )

        if B != 0:
            print(
                f"    odd(B)={odd_part(B)}"
            )

        print(
            f"    v2(A+B)="
            f"{vstr(A+B)} "
            f"v2(A-B)="
            f"{vstr(A-B)}"
        )

        print(
            f"    v2(gcd(A,B))="
            f"{gcd_val}"
        )

        if alpha < INF and beta < INF and alpha == beta:

            a = A >> alpha
            b = B >> beta

            print(
                f"    normalized odd pair="
                f"({a},{b})"
            )

            print(
                f"    odd pair mod 4="
                f"({a % 4},{b % 4})"
            )

            print(
                f"    odd pair mod 8="
                f"({a % 8},{b % 8})"
            )

        print(
            f"    depth="
            f"{min(alpha,beta) + int(equal)}"
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
Experiment 617 established:

    depth =
        min(v2(A),v2(B))
        + [v2(A)=v2(B)].

Experiment 618 attacks the remaining equality bit.

For:

    alpha = v2(A)
    beta  = v2(B),

the standard 2-adic identities imply:

    alpha != beta
        =>
    v2(A+B)=v2(A-B)=min(alpha,beta).

For:

    alpha = beta = t,

write:

    A=2^t a
    B=2^t b

with a,b odd.

Then:

    A+B = 2^t(a+b)
    A-B = 2^t(a-b).

Because a,b are both odd:

    a+b is even
    a-b is even.

Therefore:

    v2(A+B)>t
    v2(A-B)>t.

So equality of the two valuations is equivalent to:

    v2(A+B) > min(v2(A),v2(B))

and also equivalent to:

    v2(A-B) > min(v2(A),v2(B)).

This means the equality flag may be expressible entirely
through the sum/difference residuals.

The strongest possible collapse would therefore be:

    depth =
        min(v2(A),v2(B))
        +
        [
            v2(A+B)
            >
            min(v2(A),v2(B))
        ].

That would avoid explicitly comparing two valuations.

For the actual factor frames:

    FRAME A:
        A=p-3
        B=q+3

    FRAME B:
        A=p+1
        B=q-3.

The next experiment should therefore determine whether the
equality branch can be rewritten as a single 2-adic valuation
of a simple polynomial in p and q, rather than two residual
valuations plus an equality flag.
"""
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 618 START")
    print("=" * 90)
    print()
    print("EQUALITY-BRANCH COLLAPSE")
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

    test_equality_flag_from_basic_data(
        states
    )

    test_sum_difference_structure(
        states
    )

    test_odd_part_congruences(
        states
    )

    test_small_modulus_state(
        states
    )

    test_equality_vs_XY(
        states
    )

    test_equality_vs_prime_residues(
        states
    )

    test_gcd_structure(
        states
    )

    test_direct_equality_formula(
        states
    )

    print_examples(
        states
    )

    symbolic_summary()

    print("=" * 90)
    print("EXPERIMENT 618 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
