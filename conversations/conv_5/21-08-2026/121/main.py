#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 683
SHIFTED-TARGET / ODD-CORE DESCENT

Research direction:

    n = p*q

    FRAME A:
        c = 9
        C = q + 3

    FRAME B:
        c = 3
        C = p + 1

    H = gcd(C, n+c)

    depth = v2(gcd(2C, n+c))
          = min(v2(C)+1, v2(n+c))

New hypothesis:

    M0 = n + c
    M0 = 2^w * R0
    R0 odd

Instead of treating v2(n+c) only as a valuation,
factor the smaller odd target R0.

Then recursively repeat:

    Ri + c' = 2^wi * R(i+1)

for c' in {3, 9}.

Questions:

    1. Does the odd core R0 contain information about C?
    2. Does factorization of R0 recover odd(C)?
    3. Does gcd(C, R0) reproduce odd(H)?
    4. Does recursive shifted-target descent eventually expose
       additional information about C?
    5. Can a factor-dependent C be replaced by an iterated
       n-only descent?
    6. How much does the target shrink after each stripping step?

All arithmetic is exact integer arithmetic.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import gcd, isqrt


# ============================================================================
# CONFIGURATION
# ============================================================================

PRIME_LIMIT = 6000
MAX_RECURSION = 12
EXAMPLE_COUNT = 16
FACTOR_SAMPLE_COUNT = 20

# A/B examples from the previous experiments.
EXAMPLES_N = [
    9,
    15,
    21,
    33,
    39,
    57,
    69,
    77,
    87,
    93,
    111,
    141,
    183,
    213,
    485879,
    5579767,
]


# ============================================================================
# BASIC ARITHMETIC
# ============================================================================

def sieve_primes(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, flag in enumerate(sieve) if flag]


def v2(x: int) -> int:
    if x == 0:
        return 10**9

    x = abs(x)
    return (x & -x).bit_length() - 1


def odd_core(x: int) -> tuple[int, int]:
    """
    Return:

        x = 2^w * r

    with r odd.
    """
    if x == 0:
        return 0, 10**9

    w = v2(x)
    return x >> w, w


def odd_part(x: int) -> int:
    r, _ = odd_core(x)
    return r


# ============================================================================
# FACTORIZATION
# ============================================================================

def factor_integer(n: int, primes: list[int]) -> tuple[list[tuple[int, int]], int]:
    """
    Trial-factor n using the generated prime list.

    Returns:
        factors = [(prime, exponent), ...]
        remainder

    If remainder == 1, factorization is complete.
    """
    if n < 0:
        n = -n

    if n in (0, 1):
        return [], n

    x = n
    factors: list[tuple[int, int]] = []

    for p in primes:
        if p * p > x:
            break

        if x % p != 0:
            continue

        e = 0
        while x % p == 0:
            x //= p
            e += 1

        factors.append((p, e))

        if x == 1:
            break

    if x > 1:
        factors.append((x, 1))
        return factors, 1

    return factors, 1


def factor_to_string(factors: list[tuple[int, int]]) -> str:
    if not factors:
        return "1"

    parts = []
    for p, e in factors:
        if e == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{e}")

    return " * ".join(parts)


def factor_product(factors: list[tuple[int, int]]) -> int:
    out = 1
    for p, e in factors:
        out *= p ** e
    return out


# ============================================================================
# PRIME / SEMIPRIME GENERATION
# ============================================================================

def build_semiprimes(primes: list[int]) -> list[tuple[int, int, int]]:
    """
    Return all p <= q odd-prime products.
    """
    odd_primes = [p for p in primes if p & 1]

    states = []

    for i, p in enumerate(odd_primes):
        for q in odd_primes[i:]:
            states.append((p * q, p, q))

    return states


# ============================================================================
# FRAME MODEL
# ============================================================================

def frame_data(n: int, p: int, q: int) -> tuple[str, int, int]:
    """
    Return:

        frame
        C
        c

    Canonical frame from the previous experiments:

        Frame A:
            C = q + 3
            c = 9

        Frame B:
            C = p + 1
            c = 3

    Classification:
        n % 4 == 3 -> A
        n % 4 == 1 -> B
    """
    if n % 4 == 3:
        return "A", q + 3, 9

    return "B", p + 1, 3


def canonical_depth(C: int, shifted: int) -> int:
    return v2(gcd(2 * C, shifted))


# ============================================================================
# SHIFTED-TARGET DESCENT
# ============================================================================

def descend_once(value: int, c: int) -> tuple[int, int, int]:
    """
    Compute:

        value + c = 2^w * next_value

    Return:

        shifted,
        next_value,
        w
    """
    shifted = value + c
    next_value, w = odd_core(shifted)
    return shifted, next_value, w


def build_descent_chain(
    start: int,
    c_sequence: list[int],
) -> list[dict]:
    """
    Build a recursive odd-core chain.

    At each level:

        current + c = shifted
        shifted = 2^w * next

    The next state is the odd core.

    The c value may vary at every level.
    """
    chain: list[dict] = []

    current = start

    for level, c in enumerate(c_sequence):
        shifted, next_value, w = descend_once(current, c)

        chain.append(
            {
                "level": level,
                "current": current,
                "c": c,
                "shifted": shifted,
                "w": w,
                "next": next_value,
            }
        )

        if next_value == 1:
            break

        current = next_value

    return chain


# ============================================================================
# TEST HELPERS
# ============================================================================

def register_failure(
    failures: int,
    message: str,
    limit: int = 20,
) -> int:
    if failures < limit:
        print("    " + message)
    return failures + 1


def summarize_factors(factors: list[tuple[int, int]]) -> tuple[int, int]:
    """
    Return:

        number of distinct prime factors,
        total prime multiplicity
    """
    return len(factors), sum(e for _, e in factors)


# ============================================================================
# TEST 0
# ============================================================================

def test_baseline(states: list[tuple[int, int, int]]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE")
    print("=" * 90)

    failures = 0

    for n, p, q in states:
        frame, C, c = frame_data(n, p, q)

        shifted = n + c

        actual = canonical_depth(C, shifted)
        expected = min(v2(C) + 1, v2(shifted))

        if actual != expected:
            failures = register_failure(
                failures,
                (
                    f"n={n} frame={frame} "
                    f"actual={actual} expected={expected}"
                ),
            )

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ============================================================================
# TEST 1
# ============================================================================

def test_shifted_target_normalization(
    states: list[tuple[int, int, int]],
) -> int:
    print("=" * 90)
    print("TEST 1: SHIFTED TARGET NORMALIZATION")
    print("=" * 90)

    failures = 0

    distribution = Counter()
    core_sizes = []

    for n, p, q in states:
        frame, C, c = frame_data(n, p, q)

        shifted = n + c
        R, w = odd_core(shifted)

        # Exact identity.
        if shifted != (1 << w) * R:
            failures = register_failure(
                failures,
                (
                    f"n={n} shifted={shifted} "
                    f"w={w} R={R}"
                ),
            )
            continue

        distribution[(frame, w)] += 1
        core_sizes.append((n, R))

    print("valuation distribution:")
    for key in sorted(distribution):
        print(f"    {key}: {distribution[key]}")

    if core_sizes:
        reductions = [
            n / max(1, R)
            for n, R in core_sizes
        ]

        print()
        print(f"minimum n/core reduction={min(reductions):.6f}")
        print(f"maximum n/core reduction={max(reductions):.6f}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ============================================================================
# TEST 2
# ============================================================================

def test_odd_core_factorization(
    states: list[tuple[int, int, int]],
    primes: list[int],
) -> int:
    print("=" * 90)
    print("TEST 2: FACTORIZATION OF THE ODD SHIFTED CORE")
    print("=" * 90)

    failures = 0

    complete = 0
    incomplete = 0

    largest_core = 0
    largest_core_state = None

    factor_count_distribution = Counter()

    for n, p, q in states:
        frame, C, c = frame_data(n, p, q)

        shifted = n + c
        R, w = odd_core(shifted)

        factors, remainder = factor_integer(R, primes)

        reconstructed = factor_product(factors) * remainder

        if reconstructed != R:
            failures = register_failure(
                failures,
                (
                    f"n={n} R={R} "
                    f"reconstructed={reconstructed}"
                ),
            )

        if remainder == 1:
            complete += 1
        else:
            incomplete += 1

        distinct_count, multiplicity = summarize_factors(factors)

        factor_count_distribution[
            (frame, distinct_count, multiplicity)
        ] += 1

        if R > largest_core:
            largest_core = R
            largest_core_state = (n, p, q, frame, c, w, R)

    print(f"complete factorizations={complete}")
    print(f"incomplete factorizations={incomplete}")

    print("largest odd core:")
    if largest_core_state:
        n, p, q, frame, c, w, R = largest_core_state
        print(
            f"    n={n} p={p} q={q} "
            f"frame={frame} c={c} w={w} R={R}"
        )

    print()
    print("factor-count distribution:")
    for key, count in sorted(factor_count_distribution.items()):
        print(f"    {key}: {count}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ============================================================================
# TEST 3
# ============================================================================

def test_odd_core_vs_H(
    states: list[tuple[int, int, int]],
) -> int:
    print("=" * 90)
    print("TEST 3: ODD CORE VS CANONICAL gcd H")
    print("=" * 90)

    failures = 0

    odd_h_counter = Counter()
    gcd_counter = Counter()

    for n, p, q in states:
        frame, C, c = frame_data(n, p, q)

        shifted = n + c
        R, w = odd_core(shifted)

        H = gcd(C, shifted)

        odd_H = odd_part(H)

        observed = gcd(C, R)

        if observed != odd_H:
            failures = register_failure(
                failures,
                (
                    f"n={n} frame={frame} "
                    f"C={C} shifted={shifted} "
                    f"R={R} H={H} "
                    f"gcd(C,R)={observed} odd(H)={odd_H}"
                ),
            )

        odd_h_counter[(frame, odd_H)] += 1
        gcd_counter[(frame, observed)] += 1

    print("distinct odd(H):")
    for frame in ("A", "B"):
        values = {
            odd_h
            for (f, odd_h) in odd_h_counter
            if f == frame
        }

        print(
            f"    FRAME {frame}: "
            f"{len(values)} distinct"
        )

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ============================================================================
# TEST 4
# ============================================================================

def test_original_factor_exposure(
    states: list[tuple[int, int, int]],
) -> int:
    print("=" * 90)
    print("TEST 4: DOES THE ODD CORE CONTAIN p OR q?")
    print("=" * 90)

    failures = 0

    hits_p = 0
    hits_q = 0

    # This is deliberately an observational test.
    #
    # For p,q > 3:
    #
    #   n+c ≡ c (mod p)
    #
    # so p can only divide n+c when p | c.
    #
    # The experiment nevertheless checks the actual domain.

    for n, p, q in states:
        frame, C, c = frame_data(n, p, q)

        R, _ = odd_core(n + c)

        if p > 3 and R % p == 0:
            failures = register_failure(
                failures,
                f"unexpected p-divisibility: n={n} p={p} R={R}",
            )
        elif R % p == 0:
            hits_p += 1

        if q > 3 and R % q == 0:
            failures = register_failure(
                failures,
                f"unexpected q-divisibility: n={n} q={q} R={R}",
            )
        elif R % q == 0:
            hits_q += 1

    print(f"p divides odd core in exceptional cases={hits_p}")
    print(f"q divides odd core in exceptional cases={hits_q}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ============================================================================
# TEST 5
# ============================================================================

def test_core_gcd_information(
    states: list[tuple[int, int, int]],
) -> int:
    print("=" * 90)
    print("TEST 5: ODD CORE INFORMATION CONTENT")
    print("=" * 90)

    failures = 0

    relation_distribution = Counter()
    examples = []

    for n, p, q in states:
        frame, C, c = frame_data(n, p, q)

        shifted = n + c
        R, w = odd_core(shifted)

        H = gcd(C, shifted)

        odd_C = odd_part(C)
        odd_H = odd_part(H)

        g = gcd(R, odd_C)

        relation_distribution[
            (frame, v2(C), w, odd_H == g)
        ] += 1

        if g > 1 and len(examples) < 20:
            examples.append(
                (
                    n,
                    p,
                    q,
                    frame,
                    C,
                    shifted,
                    R,
                    odd_C,
                    odd_H,
                    g,
                )
            )

        if g != odd_H:
            failures = register_failure(
                failures,
                (
                    f"n={n} frame={frame} "
                    f"odd_C={odd_C} odd_H={odd_H} "
                    f"gcd(R,odd_C)={g}"
                ),
            )

    print("relation distribution:")
    for key, count in sorted(relation_distribution.items()):
        print(f"    {key}: {count}")

    print()
    print("examples where odd core contains odd canonical information:")
    for item in examples:
        (
            n,
            p,
            q,
            frame,
            C,
            shifted,
            R,
            odd_C,
            odd_H,
            g,
        ) = item

        print(
            f"    n={n} p={p} q={q} frame={frame}"
        )
        print(f"        C={C}")
        print(f"        n+c={shifted}")
        print(f"        odd-core={R}")
        print(f"        odd(C)={odd_C}")
        print(f"        odd(H)={odd_H}")
        print(f"        gcd(odd-core,odd(C))={g}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ============================================================================
# TEST 6
# ============================================================================

def test_recursive_descent(
    states: list[tuple[int, int, int]],
) -> int:
    print("=" * 90)
    print("TEST 6: RECURSIVE SHIFTED-TARGET DESCENT")
    print("=" * 90)

    failures = 0

    # Try both possible shifts after every odd-core reduction.
    #
    # This is intentionally exploratory:
    #
    #     R -> odd_core(R+3)
    #     R -> odd_core(R+9)
    #
    # We do NOT assume the next level belongs to the original frame.

    depth_distribution = Counter()
    nontrivial_gcd_distribution = Counter()
    first_hit_examples = []

    for n, p, q in states:
        frame, C, c0 = frame_data(n, p, q)

        current = n

        first_hit = None

        for level in range(MAX_RECURSION):
            results = []

            for c in (3, 9):
                shifted = current + c
                R, w = odd_core(shifted)

                g = gcd(C, R)

                results.append(
                    {
                        "c": c,
                        "shifted": shifted,
                        "R": R,
                        "w": w,
                        "g": g,
                    }
                )

                if g > 1 and g != 1:
                    if first_hit is None:
                        first_hit = (
                            level,
                            c,
                            current,
                            shifted,
                            R,
                            g,
                        )

            # Follow the original frame's shift on the first branch.
            shifted, next_value, w = descend_once(current, c0)

            depth_distribution[level] += 1

            if next_value == 1:
                break

            current = next_value

        if first_hit is not None and len(first_hit_examples) < 20:
            first_hit_examples.append(
                (
                    n,
                    p,
                    q,
                    frame,
                    C,
                    first_hit,
                )
            )

    print("first nontrivial gcd hits:")
    for item in first_hit_examples:
        n, p, q, frame, C, hit = item
        level, c, current, shifted, R, g = hit

        print(
            f"    n={n} p={p} q={q} frame={frame}"
        )
        print(
            f"        C={C}"
        )
        print(
            f"        level={level} "
            f"current={current} "
            f"c={c}"
        )
        print(
            f"        shifted={shifted} "
            f"odd-core={R} "
            f"gcd(C,odd-core)={g}"
        )

    print()
    print("descent levels observed:")
    for level, count in sorted(depth_distribution.items()):
        print(f"    level={level}: states={count}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ============================================================================
# TEST 7
# ============================================================================

def test_recursive_n_only_signatures(
    states: list[tuple[int, int, int]],
) -> int:
    print("=" * 90)
    print("TEST 7: RECURSIVE N-ONLY SIGNATURE COLLISIONS")
    print("=" * 90)

    failures = 0

    buckets: dict[tuple, set[int]] = defaultdict(set)

    for n, p, q in states:
        frame, C, c0 = frame_data(n, p, q)

        current = n
        signature = []

        for level in range(MAX_RECURSION):
            shifted = current + c0
            R, w = odd_core(shifted)

            signature.append((w, R & 0xFFFF))

            current = R

            if current == 1:
                break

        buckets[(frame, tuple(signature))].add(
            canonical_depth(C, n + c0)
        )

    ambiguous = 0
    examples = 0

    for key, depths in buckets.items():
        if len(depths) > 1:
            ambiguous += 1

            if examples < 20:
                print(
                    f"    signature={key} depths={sorted(depths)}"
                )
                examples += 1

    print()
    print(f"signatures={len(buckets)}")
    print(f"ambiguous={ambiguous}")
    print(f"failures={failures}")
    print()

    return failures


# ============================================================================
# EXAMPLES
# ============================================================================

def show_examples(
    states: list[tuple[int, int, int]],
    primes: list[int],
) -> None:
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = {
        n: (p, q)
        for n, p, q in states
        if n in EXAMPLES_N
    }

    shown = 0

    for n in EXAMPLES_N:
        if n not in wanted:
            continue

        p, q = wanted[n]

        frame, C, c = frame_data(n, p, q)

        shifted = n + c
        R, w = odd_core(shifted)

        H = gcd(C, shifted)
        depth = canonical_depth(C, shifted)

        factors, remainder = factor_integer(R, primes)

        print()
        print(
            f"n={n} p={p} q={q} frame={frame}"
        )
        print(f"    C={C}")
        print(f"    c={c}")
        print(f"    n+c={shifted}")
        print(f"    v2(n+c)={w}")
        print(f"    odd-core R={R}")
        print(f"    factor(R)={factor_to_string(factors)}")

        if remainder != 1:
            print(f"    factor remainder={remainder}")

        print(f"    H=gcd(C,n+c)={H}")
        print(f"    odd(H)={odd_part(H)}")
        print(f"    v2(C)={v2(C)}")
        print(f"    depth={depth}")

        print("    recursive descent:")

        current = n

        for level in range(min(MAX_RECURSION, 6)):
            shifted2, next_value, w2 = descend_once(
                current,
                c,
            )

            g = gcd(C, next_value)

            print(
                f"        L{level}: "
                f"x={current} -> "
                f"x+{c}={shifted2} = "
                f"2^{w2}*{next_value} "
                f"gcd(C,odd-core)={g}"
            )

            current = next_value

            if current == 1:
                break

        shown += 1

        if shown >= EXAMPLE_COUNT:
            break

    print()
    print(f"examples shown={shown}")
    print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    primes = sieve_primes(PRIME_LIMIT)
    states = build_semiprimes(primes)

    odd_primes = [p for p in primes if p & 1]

    print("=" * 90)
    print("EXPERIMENT 683 START")
    print("=" * 90)
    print()
    print(f"prime limit={PRIME_LIMIT}")
    print(f"odd primes={len(odd_primes)}")
    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_shifted_target_normalization(states)
    total_failures += test_odd_core_factorization(states, primes)
    total_failures += test_odd_core_vs_H(states)
    total_failures += test_original_factor_exposure(states)
    total_failures += test_core_gcd_information(states)
    total_failures += test_recursive_descent(states)
    total_failures += test_recursive_n_only_signatures(states)

    show_examples(states, primes)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()
    print("For the canonical frame shift:")
    print()
    print("    M = n + c")
    print("    M = 2^w * R")
    print("    R odd")
    print()
    print("with:")
    print()
    print("    FRAME A: c=9, C=q+3")
    print("    FRAME B: c=3, C=p+1")
    print()
    print("and:")
    print()
    print("    H = gcd(C, n+c)")
    print("    depth = v2(gcd(2C,n+c))")
    print()
    print("The new question is whether the smaller odd target")
    print()
    print("    R = (n+c) / 2^v2(n+c)")
    print()
    print("contains enough information to recover:")
    print()
    print("    odd(H)")
    print("    odd(C)")
    print("    v2(C)")
    print("    depth")
    print()
    print("The recursive experiment then applies")
    print()
    print("    x -> odd_core(x+3)")
    print("    x -> odd_core(x+9)")
    print()
    print("to determine whether repeated shifted-target descent")
    print("eventually exposes factor-dependent information.")
    print()
    print(f"TOTAL FAILURES={total_failures}")

    if total_failures == 0:
        print("STATUS=ALL EXACT ARITHMETIC TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 683 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
