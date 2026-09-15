#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 57
MIXED RESULTANT R(n,s) / FACTOR-POLYNOMIAL COUPLING
CYCLOTOMIC VS RANDOM CONTROL
NO CSV OUTPUT
==============================================================================

CORE OBJECT

    f(x) = x^2 - s*x + n
    F(x) = x^2 + x + 1

For x = p,q:

    f(x) = 0
    F(p)F(q) = resultant(f,F)

The resultant is

    R(n,s) = n^2 + n*s - n + s^2 + s + 1.

Thus for the true factorization:

    R(n,s) = F(p) * F(q).

This experiment asks whether R(n,s) provides information about the
true sum s beyond:

    D(n,s) = s^2 - 4n = (p-q)^2.

We compare:

    A. discriminant QR sieve alone
    B. mixed-resultant QR sieve alone
    C. discriminant + mixed-resultant sieve
    D. cyclotomic moduli
    E. random control moduli

The key test is CONDITIONAL INFORMATION:

    Does R(n,s) reject discriminant survivors
    more strongly than a generic random polynomial would?

No prime-pair enumeration is performed.
No CSV files are produced.
==============================================================================
"""

from __future__ import annotations

import math
import random
import time
from statistics import median, mean


SEED = 20260814

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

TARGETS = [
    (3318013, 4042603),
    (2129167, 3402323),
    (2224517, 3978749),
    (3685051, 4020281),
    (2399627, 2452649),
    (2593039, 2996527),
    (2149859, 2772097),
    (2060543, 2514401),
    (2675423, 2883973),
    (2828887, 3960137),
    (3497381, 3793241),
    (2193509, 4011353),
]

R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723
]

CONTROL = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]

MAX_OUTPUT_RESIDUES = 20


# -----------------------------------------------------------------------------
# Basic number theory
# -----------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> list[int]:
    """Return all primes p with lo <= p < hi."""
    limit = hi
    is_prime = bytearray(b"\x01") * limit
    if limit > 0:
        is_prime[0] = 0
    if limit > 1:
        is_prime[1] = 0

    root = int(math.isqrt(limit - 1))
    for p in range(2, root + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:limit:p] = b"\x00" * (
                ((limit - 1 - start) // p) + 1
            )

    return [p for p in range(lo, hi) if is_prime[p]]


def qr_table(mod: int) -> set[int]:
    """Quadratic residues modulo an odd prime, including zero."""
    return {(x * x) % mod for x in range(mod)}


def legendre_symbol(a: int, p: int) -> int:
    """Return Legendre symbol (a/p) for odd prime p."""
    a %= p
    if a == 0:
        return 0
    value = pow(a, (p - 1) // 2, p)
    return 1 if value == 1 else -1


# -----------------------------------------------------------------------------
# Algebraic objects
# -----------------------------------------------------------------------------

def discriminant(n: int, s: int) -> int:
    return s * s - 4 * n


def mixed_resultant(n: int, s: int) -> int:
    """
    Resultant of:

        x^2 - s*x + n
        x^2 + x + 1

    Explicitly:

        R(n,s) = n^2 + n*s - n + s^2 + s + 1
    """
    return n * n + n * s - n + s * s + s + 1


def F(x: int) -> int:
    return x * x + x + 1


def factor_pair_sum(p: int, q: int) -> int:
    return p + q


def verify_resultant_identity(p: int, q: int) -> bool:
    n = p * q
    s = p + q
    lhs = mixed_resultant(n, s)
    rhs = F(p) * F(q)
    return lhs == rhs


# -----------------------------------------------------------------------------
# Domain
# -----------------------------------------------------------------------------

def feasible_sum_bounds(p_lo: int, p_hi: int) -> tuple[int, int]:
    """
    Both factors are in [PRIME_LO, PRIME_HI), so:

        2*LO <= s <= 2*(HI-1)

    We restrict to even sums.
    """
    s_min = 2 * p_lo
    s_max = 2 * (p_hi - 1)
    if s_min % 2:
        s_min += 1
    if s_max % 2:
        s_max -= 1
    return s_min, s_max


def all_even_sums(s_min: int, s_max: int) -> list[int]:
    return list(range(s_min, s_max + 1, 2))


# -----------------------------------------------------------------------------
# Candidate filtering
# -----------------------------------------------------------------------------

def discriminant_filter(
    n: int,
    sums: list[int],
    moduli: list[int],
    qr: dict[int, set[int]],
) -> list[int]:
    survivors = sums

    for ell in moduli:
        table = qr[ell]
        survivors = [
            s for s in survivors
            if (s * s - 4 * n) % ell in table
        ]

    return survivors


def resultant_filter(
    n: int,
    sums: list[int],
    moduli: list[int],
    qr: dict[int, set[int]],
) -> list[int]:
    survivors = sums

    for ell in moduli:
        table = qr[ell]
        survivors = [
            s for s in survivors
            if mixed_resultant(n, s) % ell in table
        ]

    return survivors


def combined_filter(
    n: int,
    sums: list[int],
    moduli: list[int],
    qr: dict[int, set[int]],
) -> list[int]:
    survivors = sums

    for ell in moduli:
        table = qr[ell]
        next_survivors = []

        for s in survivors:
            d = (s * s - 4 * n) % ell
            if d not in table:
                continue

            r = mixed_resultant(n, s) % ell
            if r not in table:
                continue

            next_survivors.append(s)

        survivors = next_survivors

    return survivors


# -----------------------------------------------------------------------------
# More selective diagnostics
# -----------------------------------------------------------------------------

def resultant_zero_hits(
    n: int,
    sums: list[int],
    moduli: list[int],
) -> dict[int, int]:
    """
    Count candidates with R(n,s) == 0 mod ell.
    """
    out: dict[int, int] = {}
    for ell in moduli:
        out[ell] = sum(
            1 for s in sums
            if mixed_resultant(n, s) % ell == 0
        )
    return out


def resultant_qr_fraction(
    n: int,
    sums: list[int],
    ell: int,
    qr: set[int],
) -> float:
    if not sums:
        return 0.0
    return sum(
        1 for s in sums
        if mixed_resultant(n, s) % ell in qr
    ) / len(sums)


def discriminant_qr_fraction(
    n: int,
    sums: list[int],
    ell: int,
    qr: set[int],
) -> float:
    if not sums:
        return 0.0
    return sum(
        1 for s in sums
        if (s * s - 4 * n) % ell in qr
    ) / len(sums)


def joint_qr_fraction(
    n: int,
    sums: list[int],
    ell: int,
    qr: set[int],
) -> float:
    if not sums:
        return 0.0

    count = 0
    for s in sums:
        d = (s * s - 4 * n) % ell
        r = mixed_resultant(n, s) % ell
        if d in qr and r in qr:
            count += 1

    return count / len(sums)


def expected_independent_joint(
    d_frac: float,
    r_frac: float,
) -> float:
    return d_frac * r_frac


# -----------------------------------------------------------------------------
# Exact factor recovery from a candidate sum
# -----------------------------------------------------------------------------

def recover_factors(n: int, s: int) -> tuple[int, int] | None:
    d2 = discriminant(n, s)
    if d2 < 0:
        return None

    d = math.isqrt(d2)
    if d * d != d2:
        return None

    if (s - d) % 2 != 0:
        return None

    p = (s - d) // 2
    q = (s + d) // 2

    if p * q != n:
        return None

    return p, q


# -----------------------------------------------------------------------------
# Printing helpers
# -----------------------------------------------------------------------------

def print_header(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def fmt_bool(x: bool) -> str:
    return "PASS" if x else "FAIL"


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------

def main() -> None:
    random.seed(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 57")
    print("MIXED RESULTANT R(n,s) / FACTOR-POLYNOMIAL COUPLING")
    print("CYCLOTOMIC VS RANDOM CONTROL")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. Prime population
    # -------------------------------------------------------------------------

    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LO, PRIME_HI)
    prime_time = time.perf_counter() - t0

    print_header("1. PRIME POPULATION")
    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_time:.6f}s")

    # -------------------------------------------------------------------------
    # 2. Targets
    # -------------------------------------------------------------------------

    print_header("2. TARGETS")

    target_data = []

    for i, (p, q) in enumerate(TARGETS, 1):
        n = p * q
        s = p + q

        if p not in primes or q not in primes:
            raise RuntimeError(
                f"Target {i} factors are outside prime population."
            )

        target_data.append({
            "index": i,
            "p": p,
            "q": q,
            "n": n,
            "s": s,
        })

        print(
            f"target {i:2d}: "
            f"p={p} q={q} n={n} s={s}"
        )

    # -------------------------------------------------------------------------
    # 3. Algebraic identity
    # -------------------------------------------------------------------------

    print_header("3. RESULTANT IDENTITY")

    identity_failures = 0

    for t in target_data:
        p, q, n, s = t["p"], t["q"], t["n"], t["s"]

        lhs = mixed_resultant(n, s)
        rhs = F(p) * F(q)

        ok = lhs == rhs

        print(
            f"target {t['index']:2d}: "
            f"R(n,s)=F(p)F(q)={ok} "
            f"R mod 1000={lhs % 1000:03d}"
        )

        if not ok:
            identity_failures += 1

    print(f"identity failures = {identity_failures}")
    print(f"status = {fmt_bool(identity_failures == 0)}")

    # -------------------------------------------------------------------------
    # 4. Feasible sum domain
    # -------------------------------------------------------------------------

    s_min, s_max = feasible_sum_bounds(PRIME_LO, PRIME_HI)
    sums = all_even_sums(s_min, s_max)

    print_header("4. FEASIBLE SUM DOMAIN")
    print(f"s minimum       = {s_min:,}")
    print(f"s maximum       = {s_max:,}")
    print(f"candidate sums  = {len(sums):,}")
    print("representation  = s = s_min + 2k")

    # -------------------------------------------------------------------------
    # 5. QR tables
    # -------------------------------------------------------------------------

    all_moduli = CYCLOTOMIC + CONTROL

    qr = {}
    t0 = time.perf_counter()

    for ell in all_moduli:
        qr[ell] = qr_table(ell)

    qr_time = time.perf_counter() - t0

    print_header("5. QR TABLE PREPARATION")
    print(f"moduli = {len(all_moduli)}")
    print(f"time   = {qr_time:.6f}s")

    # -------------------------------------------------------------------------
    # 6. True-sum resultant identity and local residue diagnostics
    # -------------------------------------------------------------------------

    print_header("6. TRUE-SUM MIXED RESIDUE CHECK")

    for t in target_data:
        p, q, n, s = t["p"], t["q"], t["n"], t["s"]

        print(f"TARGET {t['index']:2d} true_s={s}")

        for ell in CYCLOTOMIC:
            d = discriminant(n, s) % ell
            r = mixed_resultant(n, s) % ell

            d_ch = legendre_symbol(d, ell)
            r_ch = legendre_symbol(r, ell)

            print(
                f"  ell={ell:5d} "
                f"D={d:5d} chi(D)={d_ch:+2d} "
                f"R={r:5d} chi(R)={r_ch:+2d}"
            )

    # -------------------------------------------------------------------------
    # 7. Baseline discriminant sieve
    # -------------------------------------------------------------------------

    print_header("7. BASELINE DISCRIMINANT SIEVE")

    base_results = {}

    for family_name, moduli in (
        ("cyclotomic", CYCLOTOMIC),
        ("control", CONTROL),
    ):
        print()
        print(f"{family_name.upper()}")

        for t in target_data:
            n = t["n"]
            s_true = t["s"]

            t0 = time.perf_counter()

            survivors = discriminant_filter(
                n, sums, moduli, qr
            )

            elapsed = time.perf_counter() - t0
            true_survives = s_true in survivors

            base_results[(family_name, t["index"])] = survivors

            print(
                f"target {t['index']:2d}: "
                f"survivors={len(survivors):5d} "
                f"time={elapsed:.6f}s "
                f"true_survives={true_survives}"
            )

    # -------------------------------------------------------------------------
    # 8. Resultant-only sieve
    # -------------------------------------------------------------------------

    print_header("8. MIXED RESULTANT-ONLY SIEVE")

    resultant_results = {}

    for family_name, moduli in (
        ("cyclotomic", CYCLOTOMIC),
        ("control", CONTROL),
    ):
        print()
        print(f"{family_name.upper()}")

        for t in target_data:
            n = t["n"]
            s_true = t["s"]

            t0 = time.perf_counter()

            survivors = resultant_filter(
                n, sums, moduli, qr
            )

            elapsed = time.perf_counter() - t0
            true_survives = s_true in survivors

            resultant_results[(family_name, t["index"])] = survivors

            print(
                f"target {t['index']:2d}: "
                f"survivors={len(survivors):5d} "
                f"time={elapsed:.6f}s "
                f"true_survives={true_survives}"
            )

    # -------------------------------------------------------------------------
    # 9. Combined sieve
    # -------------------------------------------------------------------------

    print_header("9. COMBINED DISCRIMINANT + RESULTANT SIEVE")

    combined_results = {}

    for family_name, moduli in (
        ("cyclotomic", CYCLOTOMIC),
        ("control", CONTROL),
    ):
        print()
        print(f"{family_name.upper()}")

        for t in target_data:
            n = t["n"]
            s_true = t["s"]

            t0 = time.perf_counter()

            survivors = combined_filter(
                n, sums, moduli, qr
            )

            elapsed = time.perf_counter() - t0
            true_survives = s_true in survivors

            combined_results[(family_name, t["index"])] = survivors

            print(
                f"target {t['index']:2d}: "
                f"survivors={len(survivors):5d} "
                f"time={elapsed:.6f}s "
                f"true_survives={true_survives}"
            )

            recovered = [
                recover_factors(n, s)
                for s in survivors
            ]
            recovered = [
                x for x in recovered
                if x is not None
            ]

            correct = any(
                set(pair) == {t["p"], t["q"]}
                for pair in recovered
            )

            print(
                f"             exact_recoveries={len(recovered):3d} "
                f"correct={correct}"
            )

    # -------------------------------------------------------------------------
    # 10. Conditional information test
    # -------------------------------------------------------------------------

    print_header("10. CONDITIONAL RESULTANT INFORMATION")

    print(
        "For each cyclotomic modulus we compare:"
        "\n"
        "  D-survivor fraction"
        "\n"
        "  R-survivor fraction among D-survivors"
        "\n"
        "  expected independent fraction"
        "\n"
        "  observed / expected"
    )

    ratios_c = []
    ratios_r = []

    for t in target_data:
        n = t["n"]
        base_c = base_results[("cyclotomic", t["index"])]
        base_r = base_results[("control", t["index"])]

        print()
        print(f"TARGET {t['index']:2d}")

        for label, base_survivors, ratios in (
            ("C", base_c, ratios_c),
            ("R", base_r, ratios_r),
        ):
            print(f"  {label}:")

            for ell in (
                CYCLOTOMIC if label == "C" else CONTROL
            ):
                d_frac = discriminant_qr_fraction(
                    n, sums, ell, qr[ell]
                )

                r_frac = resultant_qr_fraction(
                    n, sums, ell, qr[ell]
                )

                joint_frac = joint_qr_fraction(
                    n, sums, ell, qr[ell]
                )

                expected = expected_independent_joint(
                    d_frac,
                    r_frac,
                )

                if expected == 0.0:
                    ratio = float("nan")
                else:
                    ratio = joint_frac / expected

                cond_r = (
                    sum(
                        1 for s in base_survivors
                        if mixed_resultant(n, s) % ell in qr[ell]
                    )
                    / len(base_survivors)
                    if base_survivors
                    else 0.0
                )

                ratios.append(ratio)

                print(
                    f"    ell={ell:5d} "
                    f"D_QR={d_frac:.4f} "
                    f"R_QR={r_frac:.4f} "
                    f"joint={joint_frac:.4f} "
                    f"cond_R={cond_r:.4f} "
                    f"ratio={ratio:.4f}"
                )

    valid_c = [x for x in ratios_c if math.isfinite(x)]
    valid_r = [x for x in ratios_r if math.isfinite(x)]

    print()
    print("conditional independence summary:")
    print(
        f"  cyclotomic median ratio = "
        f"{median(valid_c):.6f}"
    )
    print(
        f"  cyclotomic mean ratio   = "
        f"{mean(valid_c):.6f}"
    )
    print(
        f"  control median ratio    = "
        f"{median(valid_r):.6f}"
    )
    print(
        f"  control mean ratio      = "
        f"{mean(valid_r):.6f}"
    )
    print("  independence benchmark = 1.000000")

    # -------------------------------------------------------------------------
    # 11. Resultant zero structure
    # -------------------------------------------------------------------------

    print_header("11. RESULTANT ZERO STRUCTURE")

    for family_name, moduli in (
        ("cyclotomic", CYCLOTOMIC),
        ("control", CONTROL),
    ):
        print()
        print(f"{family_name.upper()}")

        for t in target_data:
            n = t["n"]
            true_s = t["s"]

            zero_hits = resultant_zero_hits(
                n,
                sums,
                moduli,
            )

            true_zero = [
                ell for ell in moduli
                if mixed_resultant(n, true_s) % ell == 0
            ]

            total_zero = sum(zero_hits.values())

            print(
                f"target {t['index']:2d}: "
                f"true_zero_moduli={true_zero} "
                f"total_candidate_zero_hits={total_zero}"
            )

    # -------------------------------------------------------------------------
    # 12. True sum rank among candidate residues
    # -------------------------------------------------------------------------

    print_header("12. TRUE-SUM RANK UNDER RESULTANT RESIDUES")

    for family_name, moduli in (
        ("cyclotomic", CYCLOTOMIC),
        ("control", CONTROL),
    ):
        print()
        print(f"{family_name.upper()}")

        ranks = []

        for t in target_data:
            n = t["n"]
            true_s = t["s"]

            # Score each candidate by how many moduli classify R as QR.
            scores = []

            for s in sums:
                score = sum(
                    1
                    for ell in moduli
                    if mixed_resultant(n, s) % ell in qr[ell]
                )
                scores.append((score, s))

            scores.sort(reverse=True)

            true_score = next(
                score for score, s in scores
                if s == true_s
            )

            rank = 1 + sum(
                1 for score, _ in scores
                if score > true_score
            )

            ranks.append(rank)

            print(
                f"target {t['index']:2d}: "
                f"true_score={true_score:2d}/{len(moduli)} "
                f"rank={rank}"
            )

        print(
            f"{family_name} median true-sum rank = "
            f"{median(ranks):.1f}"
        )
        print(
            f"{family_name} mean true-sum rank = "
            f"{mean(ranks):.2f}"
        )

    # -------------------------------------------------------------------------
    # 13. Local algebra cross-check on roots
    # -------------------------------------------------------------------------

    print_header("13. ROOT-FACTOR CROSS-CHECK")

    failures = 0

    for ell in CYCLOTOMIC:
        roots = [
            w for w in range(ell)
            if (w * w + w + 1) % ell == 0
        ]

        print(
            f"ell={ell:5d} roots={roots}"
        )

        for t in target_data:
            p, q = t["p"], t["q"]
            for w in roots:
                lhs = (
                    ((p * p + p + 1) % ell)
                    * ((q * q + q + 1) % ell)
                ) % ell

                n = t["n"]
                s = t["s"]

                rhs = mixed_resultant(n, s) % ell

                if lhs != rhs:
                    failures += 1

    print(f"cross-check failures = {failures}")
    print(f"status = {fmt_bool(failures == 0)}")

    # -------------------------------------------------------------------------
    # 14. Global summary
    # -------------------------------------------------------------------------

    print_header("14. GLOBAL SUMMARY")

    base_c = [
        len(base_results[("cyclotomic", i)])
        for i in range(1, len(TARGETS) + 1)
    ]
    base_r = [
        len(base_results[("control", i)])
        for i in range(1, len(TARGETS) + 1)
    ]

    res_c = [
        len(resultant_results[("cyclotomic", i)])
        for i in range(1, len(TARGETS) + 1)
    ]
    res_r = [
        len(resultant_results[("control", i)])
        for i in range(1, len(TARGETS) + 1)
    ]

    comb_c = [
        len(combined_results[("cyclotomic", i)])
        for i in range(1, len(TARGETS) + 1)
    ]
    comb_r = [
        len(combined_results[("control", i)])
        for i in range(1, len(TARGETS) + 1)
    ]

    print(
        f"cyclotomic discriminant median = "
        f"{median(base_c):.1f}"
    )
    print(
        f"control discriminant median     = "
        f"{median(base_r):.1f}"
    )
    print(
        f"cyclotomic resultant median     = "
        f"{median(res_c):.1f}"
    )
    print(
        f"control resultant median        = "
        f"{median(res_r):.1f}"
    )
    print(
        f"cyclotomic combined median      = "
        f"{median(comb_c):.1f}"
    )
    print(
        f"control combined median         = "
        f"{median(comb_r):.1f}"
    )

    print()
    print("relative additional reduction:")
    print(
        f"  cyclotomic: "
        f"{median(base_c) / median(comb_c):.3f}x"
    )
    print(
        f"  control:    "
        f"{median(base_r) / median(comb_r):.3f}x"
    )

    print()
    print("expected interpretation:")
    print("""
    RESULTANT IDENTICALITY:
        R(n,s) = F(p)F(q) should always hold for the true factor pair.

    POSITIVE:
        The mixed resultant rejects discriminant survivors more strongly
        than expected from independent/generic QR behaviour, especially
        for cyclotomic moduli but not controls.

    NEGATIVE:
        R(n,s) behaves like another generic quadratic-residue test.

    STRONGER POSITIVE:
        The true sum ranks unusually high under the R(n,s) signature,
        and the effect survives comparison with random control moduli.

    IMPORTANT:
        A zero resultant condition R(n,s)=0 only means
            F(p)=0 or F(q)=0 modulo ell.
        It is useful only if such events occur systematically.

    The primary statistic is therefore NOT the raw number of survivors.
    It is the conditional deviation from independence:
        P(R QR | D QR) / [P(R QR) * P(D QR)]
    together with the true-sum rank.
    """)

    print("=" * 78)
    print("EXPERIMENT 57 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

