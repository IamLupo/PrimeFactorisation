#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 576
# ==============================================================================
# NORMALIZATION AUDIT + ALGEBRAIC MODULAR CHILD-BIT RULES
#
# Purpose
# -------
#
# Experiment 575 showed:
#
#     mod-4 -> mod-8 -> mod-16 -> ...
#
# has a very clean binary branching structure.
#
# However, a deeper consistency problem appeared:
#
# for r != 1 the generated x,y coordinates do not make
#
#     P = p
#     Q = q
#
# under the formulas
#
#     P = (2y-2x+3)/r
#     Q = (2y+2x-3)/r.
#
# Therefore we must stop treating the previous P,Q/n relation as established.
#
#
# This experiment performs a normalization audit first.
#
#
# Original coordinate expressions
# --------------------------------
#
#     A = 2y - 2x + 3
#     B = 2y + 2x - 3
#
# and:
#
#     P = A/r
#     Q = B/r.
#
#
# If P=p and Q=q, then:
#
#     A = r*p
#     B = r*q.
#
# Solving gives:
#
#     y = r(p+q)/4
#
#     x = (r(q-p)+6)/4.
#
#
# Experiment 575 instead used:
#
#     y = (p+q)/(4r)
#     x = (q-p+6)/(4r),
#
# which is a different normalization.
#
# We therefore test both constructions and determine exactly which one
# is algebraically compatible with P=p, Q=q.
#
#
# Main modular question
# ---------------------
#
# Once the valid normalization is isolated, test:
#
#     r_k = n mod 2^k
#
#     b_k = floor(n / 2^k) mod 2
#
# and search for exact parent-conditioned rules:
#
#     b_k = f(x mod 2^a, y mod 2^b, r)
#
# but only after verifying the underlying n/x/y relation.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 500_000

MIN_LEVEL = 2
MAX_LEVEL = 8

SHOW_RULES = True
MAX_RULE_STATES = 32


# ==============================================================================
# SAMPLE
# ==============================================================================

@dataclass(frozen=True)
class Sample:

    n: int
    p: int
    q: int

    r: int

    x_old: int
    y_old: int

    x_corr: int
    y_corr: int

    P_old: Fraction
    Q_old: Fraction

    P_corr: Fraction
    Q_corr: Fraction


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int):

    sieve = bytearray(
        b"\x01"
    ) * (
        limit + 1
    )

    if limit >= 0:
        sieve[0] = 0

    if limit >= 1:
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
# PREVIOUS NORMALIZATION
# ==============================================================================

def old_xy(
    p,
    q,
    r,
):

    denominator = 4 * r

    x_num = q - p + 6
    y_num = p + q

    if x_num % denominator != 0:
        return None

    if y_num % denominator != 0:
        return None

    x = x_num // denominator
    y = y_num // denominator

    return x, y


# ==============================================================================
# CORRECT NORMALIZATION FOR P=p, Q=q
# ==============================================================================

def corrected_xy(
    p,
    q,
    r,
):

    #
    # Require:
    #
    #     2y - 2x + 3 = r*p
    #
    #     2y + 2x - 3 = r*q
    #
    # Adding:
    #
    #     4y = r(p+q)
    #
    # Subtracting:
    #
    #     4x = r(q-p)+6
    #

    y_num = r * (
        p + q
    )

    x_num = (
        r * (q - p)
        + 6
    )

    if y_num % 4 != 0:
        return None

    if x_num % 4 != 0:
        return None

    y = y_num // 4
    x = x_num // 4

    return x, y


# ==============================================================================
# P,Q
# ==============================================================================

def pq_from_xy(
    x,
    y,
    r,
):

    A = (
        2 * y
        - 2 * x
        + 3
    )

    B = (
        2 * y
        + 2 * x
        - 3
    )

    P = Fraction(
        A,
        r,
    )

    Q = Fraction(
        B,
        r,
    )

    return P, Q


# ==============================================================================
# GENERATE SAMPLES
# ==============================================================================

def generate_samples():

    sieve = prime_sieve(
        MAX_N
    )

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    samples = []

    for n, p, q in semiprimes:

        for r in (
            1,
            2,
            4,
        ):

            old = old_xy(
                p,
                q,
                r,
            )

            corr = corrected_xy(
                p,
                q,
                r,
            )

            if old is None and corr is None:
                continue

            if old is None:
                ox = oy = None
                P_old = Q_old = None

            else:
                ox, oy = old

                P_old, Q_old = pq_from_xy(
                    ox,
                    oy,
                    r,
                )

            if corr is None:
                cx = cy = None
                P_corr = Q_corr = None

            else:

                cx, cy = corr

                P_corr, Q_corr = pq_from_xy(
                    cx,
                    cy,
                    r,
                )

            samples.append(
                Sample(
                    n=n,
                    p=p,
                    q=q,
                    r=r,
                    x_old=ox,
                    y_old=oy,
                    x_corr=cx,
                    y_corr=cy,
                    P_old=P_old,
                    Q_old=Q_old,
                    P_corr=P_corr,
                    Q_corr=Q_corr,
                )
            )

    return samples


# ==============================================================================
# NORMALIZATION AUDIT
# ==============================================================================

def normalization_audit(
    samples,
):

    print()
    print("=" * 90)
    print(
        "NORMALIZATION AUDIT"
    )
    print("=" * 90)

    old_total = 0
    old_exact = 0

    corr_total = 0
    corr_exact = 0

    by_r = defaultdict(
        lambda: {
            "old_total": 0,
            "old_exact": 0,
            "corr_total": 0,
            "corr_exact": 0,
        }
    )

    examples = []

    for s in samples:

        rinfo = by_r[s.r]

        if s.x_old is not None:

            old_total += 1
            rinfo["old_total"] += 1

            if (
                s.P_old == s.p
                and s.Q_old == s.q
            ):

                old_exact += 1
                rinfo["old_exact"] += 1

            elif len(examples) < 8:

                examples.append(
                    (
                        "OLD",
                        s,
                    )
                )

        if s.x_corr is not None:

            corr_total += 1
            rinfo["corr_total"] += 1

            if (
                s.P_corr == s.p
                and s.Q_corr == s.q
            ):

                corr_exact += 1
                rinfo["corr_exact"] += 1

    print()
    print(
        "OLD NORMALIZATION"
    )

    print(
        f"    exact P=p,Q=q: "
        f"{old_exact}/{old_total}"
    )

    print()
    print(
        "CORRECTED NORMALIZATION"
    )

    print(
        f"    exact P=p,Q=q: "
        f"{corr_exact}/{corr_total}"
    )

    print()
    print(
        "BY r"
    )

    for r in sorted(by_r):

        info = by_r[r]

        print(
            f"    r={r}"
        )

        print(
            f"        old: "
            f"{info['old_exact']}/"
            f"{info['old_total']}"
        )

        print(
            f"        corr: "
            f"{info['corr_exact']}/"
            f"{info['corr_total']}"
        )

    print()
    print(
        "COUNTEREXAMPLES TO OLD NORMALIZATION"
    )

    for kind, s in examples:

        print(
            f"    n={s.n} "
            f"p={s.p} "
            f"q={s.q} "
            f"r={s.r} "
            f"x_old={s.x_old} "
            f"y_old={s.y_old} "
            f"P_old={s.P_old} "
            f"Q_old={s.Q_old}"
        )


# ==============================================================================
# CORRECTED SAMPLE FILTER
# ==============================================================================

def corrected_samples(
    samples,
):

    result = []

    for s in samples:

        if s.x_corr is None:
            continue

        if (
            s.P_corr != s.p
            or s.Q_corr != s.q
        ):
            continue

        result.append(
            s
        )

    return result


# ==============================================================================
# MODULAR HELPERS
# ==============================================================================

def modulus(level):
    return 1 << level


def parent_residue(
    n,
    level,
):
    return n % modulus(level)


def next_bit(
    n,
    level,
):
    return (
        n >> level
    ) & 1


# ==============================================================================
# GROUPING
# ==============================================================================

def group_parent(
    samples,
    level,
):

    groups = defaultdict(list)

    for s in samples:

        groups[
            parent_residue(
                s.n,
                level,
            )
        ].append(s)

    return groups


# ==============================================================================
# FEATURE FUNCTIONS
# ==============================================================================

def corrected_x(
    sample,
    m,
):
    return sample.x_corr % m


def corrected_y(
    sample,
    m,
):
    return sample.y_corr % m


def p_mod(
    sample,
    m,
):
    return sample.P_corr % m


def q_mod(
    sample,
    m,
):
    return sample.Q_corr % m


# ==============================================================================
# ALGEBRAIC n FROM CORRECTED x,y,r
# ==============================================================================

def reconstructed_n(
    sample,
):

    x = sample.x_corr
    y = sample.y_corr
    r = sample.r

    #
    # From:
    #
    #     r*p = 2y-2x+3
    #     r*q = 2y+2x-3
    #
    # multiply:
    #
    #     r^2*p*q
    #
    # = (2y-2x+3)(2y+2x-3)
    #
    # = 4y^2 - (2x-3)^2
    #
    #

    numerator = (
        4 * y * y
        - (2 * x - 3) ** 2
    )

    denominator = r * r

    if numerator % denominator != 0:
        return None

    return numerator // denominator


# ==============================================================================
# N RECONSTRUCTION AUDIT
# ==============================================================================

def reconstruction_audit(
    samples,
):

    print()
    print("=" * 90)
    print(
        "n RECONSTRUCTION AUDIT"
    )
    print("=" * 90)

    total = 0
    exact = 0

    by_r = Counter()

    for s in samples:

        n2 = reconstructed_n(
            s
        )

        if n2 is None:
            continue

        total += 1
        by_r[s.r] += 1

        if n2 == s.n:
            exact += 1

    print(
        f"n reconstructed exactly: "
        f"{exact}/{total}"
    )

    for r in sorted(by_r):

        print(
            f"    r={r}: "
            f"{by_r[r]}"
        )


# ==============================================================================
# EXACT CHILD-BIT RULE
# ==============================================================================

def exact_rule(
    subset,
    feature_function,
    feature_modulus,
    level,
):

    table = defaultdict(set)

    for s in subset:

        key = feature_function(
            s,
            feature_modulus,
        )

        bit = next_bit(
            s.n,
            level,
        )

        table[key].add(
            bit
        )

    if not table:
        return None

    for values in table.values():

        if len(values) != 1:
            return None

    return {
        key: next(iter(values))
        for key, values in table.items()
    }


# ==============================================================================
# BRANCH STRUCTURE
# ==============================================================================

def branch_structure(
    samples,
):

    print()
    print("=" * 90)
    print(
        "CORRECTED BRANCH STRUCTURE"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        groups = group_parent(
            samples,
            level,
        )

        both = 0

        print()
        print(
            f"MOD {modulus(level)} "
            f"-> MOD {modulus(level + 1)}"
        )

        for parent in sorted(groups):

            subset = groups[parent]

            bits = {
                next_bit(
                    s.n,
                    level,
                )
                for s in subset
            }

            if bits == {0, 1}:
                both += 1

            print(
                f"    parent={parent:<4} "
                f"count={len(subset):<7} "
                f"children={sorted(bits)}"
            )

        print(
            "    both-child parents =",
            both,
        )


# ==============================================================================
# DIRECT ALGEBRAIC BIT TEST
# ==============================================================================

def direct_algebraic_bit(
    sample,
    level,
):

    """
    Since reconstructed_n(sample) == n,
    compute the child bit from the corrected x,y,r directly.

        n = [4y^2-(2x-3)^2]/r^2

    """

    n2 = reconstructed_n(
        sample
    )

    if n2 is None:
        return None

    return (
        n2 >> level
    ) & 1


def algebraic_bit_audit(
    samples,
):

    print()
    print("=" * 90)
    print(
        "DIRECT ALGEBRAIC CHILD-BIT AUDIT"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        total = 0
        exact = 0

        for s in samples:

            b1 = next_bit(
                s.n,
                level,
            )

            b2 = direct_algebraic_bit(
                s,
                level,
            )

            if b2 is None:
                continue

            total += 1

            if b1 == b2:
                exact += 1

        print(
            f"MOD {modulus(level)} "
            f"next bit: "
            f"{exact}/{total}"
        )


# ==============================================================================
# PARENT-CONDITIONED RULE SEARCH
# ==============================================================================

def search_coordinate_rules(
    samples,
):

    print()
    print("=" * 90)
    print(
        "PARENT-CONDITIONED COORDINATE RULES"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_modulus = modulus(
            level
        )

        feature_modulus = modulus(
            level - 1
        )

        groups = group_parent(
            samples,
            level,
        )

        print()
        print(
            f"MOD {parent_modulus} "
            f"-> MOD {parent_modulus * 2}"
        )

        print(
            f"    expected x modulus = "
            f"{feature_modulus}"
        )

        exact_x = 0
        exact_y = 0
        exact_p = 0
        exact_q = 0

        total = 0

        for parent in sorted(groups):

            subset = groups[parent]

            bits = {
                next_bit(
                    s.n,
                    level,
                )
                for s in subset
            }

            if bits != {0, 1}:
                continue

            total += 1

            rules = {
                "x": exact_rule(
                    subset,
                    corrected_x,
                    feature_modulus,
                    level,
                ),
                "y": exact_rule(
                    subset,
                    corrected_y,
                    feature_modulus,
                    level,
                ),
                "P": exact_rule(
                    subset,
                    p_mod,
                    feature_modulus,
                    level,
                ),
                "Q": exact_rule(
                    subset,
                    q_mod,
                    feature_modulus,
                    level,
                ),
            }

            if rules["x"] is not None:
                exact_x += 1

            if rules["y"] is not None:
                exact_y += 1

            if rules["P"] is not None:
                exact_p += 1

            if rules["Q"] is not None:
                exact_q += 1

            if SHOW_RULES:

                for name, rule in rules.items():

                    if rule is None:
                        continue

                    states = len(
                        rule
                    )

                    ones = [
                        key
                        for key, bit in sorted(
                            rule.items(),
                            key=lambda item: str(item[0]),
                        )
                        if bit == 1
                    ]

                    if len(ones) <= MAX_RULE_STATES:

                        print(
                            f"    parent={parent:<4} "
                            f"{name}: "
                            f"states={states:<3} "
                            f"bit1={ones}"
                        )

        print(
            f"    x exact = {exact_x}/{total}"
        )

        print(
            f"    y exact = {exact_y}/{total}"
        )

        print(
            f"    P exact = {exact_p}/{total}"
        )

        print(
            f"    Q exact = {exact_q}/{total}"
        )


# ==============================================================================
# LOW-LEVEL ALGEBRAIC REDUCTION
# ==============================================================================

def parity_state(
    sample,
    level,
):

    """
    Compute the numerator:

        M = 4y^2 - (2x-3)^2

    and inspect its residue modulo powers of two.

    Since:

        n = M / r^2,

    the parity structure of M depends strongly on r.
    """

    x = sample.x_corr
    y = sample.y_corr
    r = sample.r

    M = (
        4 * y * y
        - (2 * x - 3) ** 2
    )

    return M, r


# ==============================================================================
# RESIDUE TABLE BY r
# ==============================================================================

def algebraic_residue_table(
    samples,
):

    print()
    print("=" * 90)
    print(
        "ALGEBRAIC RESIDUE TABLE"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        m = modulus(level)

        print()
        print(
            f"MOD {m}"
        )

        by_r = defaultdict(
            Counter
        )

        for s in samples:

            M, r = parity_state(
                s,
                level,
            )

            b = next_bit(
                s.n,
                level,
            )

            key = M % (
                m * r * r
            )

            by_r[r][
                (
                    key,
                    b,
                )
            ] += 1

        for r in sorted(by_r):

            print(
                f"    r={r}"
            )

            counter = by_r[r]

            for (
                key,
                bit,
            ), count in counter.most_common(
                16
            ):

                print(
                    f"        M residue={key:<6} "
                    f"bit={bit} "
                    f"count={count}"
                )


# ==============================================================================
# EXAMPLE CORRECTED STATES
# ==============================================================================

def examples(
    samples,
):

    print()
    print("=" * 90)
    print(
        "CORRECTED COORDINATE EXAMPLES"
    )
    print("=" * 90)

    seen = set()

    for s in samples:

        if s.n in seen:
            continue

        seen.add(
            s.n
        )

        print()
        print(
            f"n={s.n} "
            f"p={s.p} "
            f"q={s.q} "
            f"r={s.r}"
        )

        print(
            f"    corrected x={s.x_corr}"
        )

        print(
            f"    corrected y={s.y_corr}"
        )

        print(
            f"    P={s.P_corr}"
        )

        print(
            f"    Q={s.Q_corr}"
        )

        print(
            f"    reconstructed n="
            f"{reconstructed_n(s)}"
        )

        if len(seen) >= 10:
            break


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 576 START"
    )
    print("=" * 90)

    print()
    print(
        "NORMALIZATION AUDIT + "
        "ALGEBRAIC MODULAR CHILD-BIT RULES"
    )

    # --------------------------------------------------------------------------
    # Generate both coordinate constructions.
    # --------------------------------------------------------------------------

    samples = generate_samples()

    print()
    print(
        "Generated raw samples =",
        len(samples),
    )

    # --------------------------------------------------------------------------
    # Audit normalization.
    # --------------------------------------------------------------------------

    normalization_audit(
        samples
    )

    # --------------------------------------------------------------------------
    # Keep only algebraically consistent corrected states.
    # --------------------------------------------------------------------------

    valid = corrected_samples(
        samples
    )

    print()
    print(
        "Algebraically valid corrected samples =",
        len(valid),
    )

    if not valid:

        print(
            "No valid corrected samples."
        )

        return

    # --------------------------------------------------------------------------
    # Audit n reconstruction.
    # --------------------------------------------------------------------------

    reconstruction_audit(
        valid
    )

    # --------------------------------------------------------------------------
    # Direct algebraic child bit.
    # --------------------------------------------------------------------------

    algebraic_bit_audit(
        valid
    )

    # --------------------------------------------------------------------------
    # Branch structure.
    # --------------------------------------------------------------------------

    branch_structure(
        valid
    )

    # --------------------------------------------------------------------------
    # Coordinate rules.
    # --------------------------------------------------------------------------

    search_coordinate_rules(
        valid
    )

    # --------------------------------------------------------------------------
    # Algebraic residues.
    # --------------------------------------------------------------------------

    algebraic_residue_table(
        valid
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    examples(
        valid
    )

    # --------------------------------------------------------------------------
    # Final.
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print(
        "FINAL INTERPRETATION"
    )
    print("=" * 90)

    print(
        """
Experiment 576 deliberately separates two questions.

1. Is the x,y,r normalization algebraically consistent with P and Q?

2. Once the normalization is correct, does the next modular bit depend
   on a parent-conditioned function of x,y,r?

The exact identity being tested is:

    r^2 n
      =
    4y^2 - (2x-3)^2

provided that:

    2y-2x+3 = r*p
    2y+2x-3 = r*q.

The modular hierarchy is then studied directly from:

    n mod 2^k
        ->
    next bit
        ->
    x,y,r residue state.

This removes the normalization ambiguity present in Experiments 571-575.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 576 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
