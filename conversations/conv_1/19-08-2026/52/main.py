#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 219 — EXACT n=pq FAMILY CROSS-CASE SCHUR / 7-ADIC AUDIT
==============================================================================

Purpose
-------

The previous experiments established a very strong theorem-like structure
for the n=6=2*3 case, but that is still a single instance.

Experiment 219 changes the target:

    compare genuinely independent n=pq instances.

This script does NOT invent formulas for q1(p,q), q3(p,q), H, B, C, K,
or the retained matrix.

Instead, enter exact outputs from the actual n=pq construction in CASES.

Each case contains:

    p, q, q1, q3

and optionally the exact structural data:

    det_H,
    rank_F,
    rank_A,
    rank_K,
    H_content,
    K_content,
    C_content.

The script then checks:

    1. projective 7-adic invariance mechanisms;
    2. the general Hensel slope formula;
    3. logarithmic principal-coordinate reconstruction;
    4. common-source scaling invariance;
    5. cross-case identities actually present in the supplied family data.

A family theorem is NEVER claimed from one case.

No other Python file is imported.
"""


from __future__ import annotations

from fractions import Fraction
import sys


# ============================================================================
# EXACT CASE DATA
# ============================================================================
#
# IMPORTANT:
#
# Replace/add entries ONLY with values produced by the actual n=pq
# construction.
#
# The first row is the already-established n=6=2*3 case.
#
# When you have the next genuine case, e.g. n=10=2*5 or n=15=3*5,
# add it here.
#
# Optional structural fields may be left as None until available.
# ============================================================================

CASES = [
    {
        "p": 2,
        "q": 3,
        "q1": 29144191,
        "q3": 24794967,

        # Optional exact structural data from the real matrix construction:
        "det_H": 2,
        "rank_F": 4,
        "rank_A": 3,
        "rank_K": 1,
        "H_content": 1,
        "C_content": 17,
        "K_content": 17,
    },

    # Example placeholder ONLY.
    #
    # DO NOT fabricate these values.
    #
    # {
    #     "p": 2,
    #     "q": 5,
    #     "q1": ...,
    #     "q3": ...,
    #     "det_H": ...,
    #     "rank_F": ...,
    #     "rank_A": ...,
    #     "rank_K": ...,
    #     "H_content": ...,
    #     "C_content": ...,
    #     "K_content": ...,
    # },
]


P7 = 7
END_E = 12
PRINCIPAL_BASE = 729


# ============================================================================
# EXACT ARITHMETIC
# ============================================================================

def valuation_p(x: int, p: int) -> int | None:
    x = int(x)
    p = int(p)

    if p < 2:
        raise ValueError("p must be >= 2.")

    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % p == 0:
        x //= p
        v += 1

    return v


def valuation_7(x: int) -> int | None:
    return valuation_p(x, 7)


def gcd_int(a: int, b: int) -> int:
    a = abs(int(a))
    b = abs(int(b))

    while b:
        a, b = b, a % b

    return a


def inverse_mod(a: int, m: int) -> int:
    a = int(a)
    m = int(m)

    if m <= 1:
        return 0

    a %= m

    if a == 0:
        raise ArithmeticError(
            f"Cannot invert 0 modulo {m}."
        )

    old_r, r = a, m
    old_s, s = 1, 0

    while r:
        q = old_r // r

        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s

    if old_r != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return int(old_s % m)


def ratio_mod(
    q1: int,
    q3: int,
    modulus: int,
) -> int:

    q1 = int(q1)
    q3 = int(q3)
    modulus = int(modulus)

    return int(
        (
            q1
            * inverse_mod(q3, modulus)
        )
        % modulus
    )


def orbit_mod(
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            2
            * pow(
                3,
                int(k),
                int(modulus),
            )
        )
        % int(modulus)
    )


# ============================================================================
# HENSEL RECONSTRUCTION
# ============================================================================

def find_k1(
    q1: int,
    q3: int,
) -> int:

    rho = ratio_mod(q1, q3, 7)

    matches = []

    for k in range(6):

        if orbit_mod(k, 7) == rho:
            matches.append(k)

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected unique k1, got {matches} for rho={rho}."
        )

    return int(matches[0])


def audit_hensel_case(
    q1: int,
    q3: int,
    end_e: int,
):

    k1 = find_k1(q1, q3)

    current_k = int(k1)

    exponents = [current_k]
    digits = []
    slopes = []
    intercepts = []

    for e in range(1, end_e):

        modulus_e = int(P7 ** e)
        modulus_next = int(P7 ** (e + 1))

        target = ratio_mod(
            q1,
            q3,
            modulus_next,
        )

        order = int(
            6 * P7 ** (e - 1)
        )

        candidates = [
            int(current_k + t * order)
            for t in range(7)
        ]

        residuals = []

        for candidate in candidates:

            orbit = orbit_mod(
                candidate,
                modulus_next,
            )

            residual = int(
                (orbit - target)
                % modulus_next
            )

            if residual % modulus_e != 0:
                raise ArithmeticError(
                    f"Candidate {candidate} did not descend "
                    f"from level e={e}."
                )

            residuals.append(
                int(
                    (residual // modulus_e)
                    % P7
                )
            )

        A = residuals[0]

        B = int(
            (residuals[1] - residuals[0])
            % P7
        )

        affine = all(
            residuals[t]
            == (A + B * t) % P7
            for t in range(7)
        )

        matches = [
            t
            for t in range(7)
            if orbit_mod(
                candidates[t],
                modulus_next,
            ) == target
        ]

        if len(matches) != 1:
            raise ArithmeticError(
                f"Expected unique lift at e={e}, "
                f"got matches={matches}."
            )

        chosen_t = int(matches[0])
        next_k = int(candidates[chosen_t])

        digits.append(chosen_t)
        slopes.append(B)
        intercepts.append(A)
        exponents.append(next_k)

        current_k = next_k

        if not affine:
            raise ArithmeticError(
                f"Non-affine residual law at e={e}."
            )

    rho7 = ratio_mod(q1, q3, 7)

    predicted_slope = int(
        (
            rho7
            * (
                (PRINCIPAL_BASE - 1)
                // 7
            )
        )
        % 7
    )

    slope_formula_exact = all(
        B == predicted_slope
        for B in slopes
    )

    return {
        "k1": k1,
        "rho7": rho7,
        "exponents": exponents,
        "digits": digits,
        "slopes": slopes,
        "intercepts": intercepts,
        "predicted_slope": predicted_slope,
        "slope_formula_exact": slope_formula_exact,
    }


# ============================================================================
# 7-ADIC LOGARITHM
# ============================================================================

def fraction_valuation_7(
    x: Fraction,
) -> int | None:

    if x == 0:
        return None

    vn = valuation_7(
        x.numerator
    )

    vd = valuation_7(
        x.denominator
    )

    vn = 0 if vn is None else vn
    vd = 0 if vd is None else vd

    return vn - vd


def log_one_plus_x(
    x: Fraction,
    precision: int,
) -> Fraction:

    vx = fraction_valuation_7(x)

    if vx is None or vx <= 0:
        raise ArithmeticError(
            "7-adic logarithm requires v7(x)>0."
        )

    total = Fraction(0)
    n = 1

    while True:

        term = (
            x ** n
        ) / n

        tv = fraction_valuation_7(
            term
        )

        sign = (
            1
            if n % 2 == 1
            else -1
        )

        total += sign * term

        if tv is not None and tv >= precision:

            safe = True

            for j in range(
                n + 1,
                n + 8,
            ):

                future = (
                    j * vx
                    - (
                        valuation_7(j)
                        or 0
                    )
                )

                if future < precision:
                    safe = False
                    break

            if safe:
                break

        n += 1

        if n > 4 * precision + 100:
            raise ArithmeticError(
                "Logarithm truncation exceeded safe bound."
            )

    return total


def fraction_mod(
    x: Fraction,
    modulus: int,
) -> int:

    modulus = int(modulus)

    if x.denominator % 7 == 0:
        raise ArithmeticError(
            "Rational denominator is divisible by 7."
        )

    return int(
        (
            (x.numerator % modulus)
            * inverse_mod(
                x.denominator % modulus,
                modulus,
            )
        )
        % modulus
    )


def logarithmic_coordinate(
    q1: int,
    q3: int,
    k1: int,
):

    base = int(
        2 * 3 ** int(k1)
    )

    rho = Fraction(
        int(q1),
        int(q3),
    )

    argument = (
        rho / base
    ) - 1

    if (
        fraction_valuation_7(
            argument
        )
        is None
        or
        fraction_valuation_7(
            argument
        ) <= 0
    ):
        return {
            "defined": False,
            "m": None,
        }

    log_source = log_one_plus_x(
        argument,
        END_E + 4,
    )

    log_base = log_one_plus_x(
        Fraction(
            PRINCIPAL_BASE - 1,
            1,
        ),
        END_E + 4,
    )

    logarithmic_m = (
        log_source
        / log_base
    )

    return {
        "defined": True,
        "m": logarithmic_m,
    }


# ============================================================================
# TERMINAL SCHUR SOURCE ROW
# ============================================================================

def terminal_row(
    q1: int,
    q3: int,
):
    return [
        int(q1 - 2 * q3),
        int(q1 - 6 * q3),
        int(q1 - 18 * q3),
    ]


def row_gcd(
    row,
) -> int:

    values = [
        abs(int(x))
        for x in row
        if int(x) != 0
    ]

    if not values:
        return 0

    g = values[0]

    for value in values[1:]:
        g = gcd_int(
            g,
            value,
        )

    return g


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 219 — EXACT n=pq FAMILY CROSS-CASE "
        "SCHUR / 7-ADIC AUDIT"
    )
    print("=" * 78)

    if not CASES:
        raise RuntimeError(
            "CASES is empty."
        )

    # ------------------------------------------------------------------
    # 1. CASE AUDIT
    # ------------------------------------------------------------------

    results = []

    print()
    print("=" * 78)
    print("1. CASE-BY-CASE EXACT AUDIT")
    print("=" * 78)

    for idx, case in enumerate(
        CASES,
        start=1,
    ):

        p = int(case["p"])
        q = int(case["q"])
        q1 = int(case["q1"])
        q3 = int(case["q3"])

        if p <= 1 or q <= 1:
            raise ValueError(
                f"Invalid prime parameters in case {idx}."
            )

        n = int(
            p * q
        )

        print()
        print(
            f"  CASE {idx}: "
            f"(p,q)=({p},{q}) "
            f"n={n}"
        )

        print(
            f"    q1={q1}"
        )

        print(
            f"    q3={q3}"
        )

        print(
            f"    gcd(q1,q3)="
            f"{gcd_int(q1,q3)}"
        )

        hensel = audit_hensel_case(
            q1,
            q3,
            END_E,
        )

        log_data = logarithmic_coordinate(
            q1,
            q3,
            hensel["k1"],
        )

        modulus = int(
            7 ** (END_E - 1)
        )

        if log_data["defined"]:

            log_m_residue = int(
                fraction_mod(
                    log_data["m"],
                    modulus,
                )
            )

            principal_m = int(
                (
                    hensel["exponents"][-1]
                    - hensel["k1"]
                )
                // 6
            )

            log_match = (
                log_m_residue
                == (
                    principal_m
                    % modulus
                )
            )

        else:

            log_m_residue = None
            principal_m = int(
                (
                    hensel["exponents"][-1]
                    - hensel["k1"]
                )
                // 6
            )
            log_match = False

        row = terminal_row(
            q1,
            q3,
        )

        structural = {
            "det_H": case.get("det_H"),
            "rank_F": case.get("rank_F"),
            "rank_A": case.get("rank_A"),
            "rank_K": case.get("rank_K"),
            "H_content": case.get("H_content"),
            "C_content": case.get("C_content"),
            "K_content": case.get("K_content"),
        }

        result = {
            "p": p,
            "q": q,
            "n": n,
            "q1": q1,
            "q3": q3,
            "rho7": hensel["rho7"],
            "k1": hensel["k1"],
            "final_k": hensel["exponents"][-1],
            "principal_m": principal_m,
            "digits": hensel["digits"],
            "slopes": hensel["slopes"],
            "predicted_slope": hensel["predicted_slope"],
            "slope_formula_exact": hensel[
                "slope_formula_exact"
            ],
            "log_defined": log_data["defined"],
            "log_m_residue": log_m_residue,
            "log_match": log_match,
            "terminal_row": row,
            "terminal_row_gcd": row_gcd(row),
            "structural": structural,
        }

        results.append(
            result
        )

        print(
            f"    rho_mod7={result['rho7']}"
        )

        print(
            f"    k1={result['k1']}"
        )

        print(
            f"    final_k={result['final_k']}"
        )

        print(
            f"    principal_m={result['principal_m']}"
        )

        print(
            f"    slope_sequence={result['slopes']}"
        )

        print(
            f"    predicted_slope="
            f"{result['predicted_slope']}"
        )

        print(
            f"    slope_formula_exact="
            f"{result['slope_formula_exact']}"
        )

        print(
            f"    logarithmic_match="
            f"{result['log_match']}"
        )

        print(
            f"    terminal_row_gcd="
            f"{result['terminal_row_gcd']}"
        )

        for name, value in structural.items():

            if value is not None:
                print(
                    f"    {name}={value}"
                )

    # ------------------------------------------------------------------
    # 2. CROSS-CASE COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. CROSS-CASE COMPARISON")
    print("=" * 78)

    if len(results) < 2:

        print(
            "  family_size=1"
        )

        print(
            "  cross_case_comparison=NOT_YET_AVAILABLE"
        )

    else:

        fields = [
            "rho7",
            "k1",
            "predicted_slope",
            "slope_formula_exact",
            "log_match",
        ]

        for field in fields:

            values = [
                result[field]
                for result in results
            ]

            print(
                f"  {field}={values}"
            )

        structural_fields = [
            "det_H",
            "rank_F",
            "rank_A",
            "rank_K",
            "H_content",
            "C_content",
            "K_content",
        ]

        for field in structural_fields:

            values = [
                result["structural"][field]
                for result in results
                if result["structural"][field] is not None
            ]

            if values:
                print(
                    f"  {field}={values}"
                )

    # ------------------------------------------------------------------
    # 3. CANDIDATE UNIVERSALITIES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CANDIDATE UNIVERSAL FAMILY IDENTITIES")
    print("=" * 78)

    all_slope_exact = all(
        result["slope_formula_exact"]
        for result in results
    )

    all_log_exact = all(
        result["log_match"]
        for result in results
    )

    known_det_H = [
        result["structural"]["det_H"]
        for result in results
        if result["structural"]["det_H"] is not None
    ]

    known_rank_drop = []

    for result in results:

        rank_f = result["structural"]["rank_F"]
        rank_a = result["structural"]["rank_A"]

        if (
            rank_f is not None
            and rank_a is not None
        ):
            known_rank_drop.append(
                int(rank_f - rank_a)
            )

    known_rank_k = [
        result["structural"]["rank_K"]
        for result in results
        if result["structural"]["rank_K"] is not None
    ]

    print(
        f"  all_slope_formulas_exact="
        f"{all_slope_exact}"
    )

    print(
        f"  all_logarithmic_matches="
        f"{all_log_exact}"
    )

    print(
        f"  known_det_H_values={known_det_H}"
    )

    print(
        f"  known_rank_drop_values="
        f"{known_rank_drop}"
    )

    print(
        f"  known_rank_K_values="
        f"{known_rank_k}"
    )

    # ------------------------------------------------------------------
    # 4. FAMILY STATUS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. FAMILY-THEOREM STATUS")
    print("=" * 78)

    family_size = len(results)

    if family_size < 2:

        print(
            "  family_generalization_established=False"
        )

        print(
            "  reason="
            "only one genuine n=pq case supplied"
        )

    else:

        print(
            "  family_generalization_candidate=True"
        )

        print(
            f"  independent_cases={family_size}"
        )

    # ------------------------------------------------------------------
    # 5. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The purpose of Experiment 219 is to cross the boundary from:

    one exactly analyzed n=pq instance

to:

    several genuinely independent n=pq instances.

The 7-adic part can already be expressed from each source ratio:

    rho = q1/q3,

    k = k1 + 6m,

    B = rho*(729-1)/7 mod 7.

The unresolved family questions are on the source side:

    q1 = q1(p,q),
    q3 = q3(p,q),

and on the Schur side:

    H = H(p,q),
    B = B(p,q),
    C = C(p,q),
    K = K(p,q).

Only independent exact n=pq cases can establish whether the observed
Schur rank drop, boundary determinant, terminal content, and source
7-adic structures are universal.

This experiment therefore treats every supplied case equally and makes
no family claim from a single row.
"""
    )

    # ------------------------------------------------------------------
    # 6. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(results) >= 1
        and all_slope_exact
        and all_log_exact
    )

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  cases_nonempty="
        f"{len(results) >= 1}"
    )

    print(
        f"  all_case_slope_formulas_exact="
        f"{all_slope_exact}"
    )

    print(
        f"  all_case_logarithmic_reconstructions_exact="
        f"{all_log_exact}"
    )

    print(
        f"  family_size={family_size}"
    )

    print(
        f"  family_theorem_yet_unproved="
        f"{family_size < 2}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 219 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise

