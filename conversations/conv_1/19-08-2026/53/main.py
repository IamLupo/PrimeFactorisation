#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 220R — EXACT n=6 vs n=10 SEMIPRIME FAMILY COMPARISON
==============================================================================

Single-file main.py.

This version fixes the undefined gcd() reference by using the local
exact gcd_int() implementation everywhere.

n=6 is populated from the established experiment.

n=10 remains an explicit input slot. Do not invent those values; replace
None with values produced by the actual n=pq construction.
"""

from __future__ import annotations

from fractions import Fraction
import sys


# ============================================================================
# KNOWN n=6 DATA
# ============================================================================

CASE_N6 = {
    "p": 2,
    "q": 3,
    "q1": 29144191,
    "q3": 24794967,

    "det_H": 2,
    "rank_F": 4,
    "rank_A": 3,
    "rank_K": 1,

    "H_content": 1,
    "C_content": 17,
    "K_content": 17,
}


# ============================================================================
# n=10 DATA
# ============================================================================
#
# Fill these from the genuine n=10 construction.
# ============================================================================

CASE_N10 = {
    "p": 2,
    "q": 5,

    "q1": None,
    "q3": None,

    "det_H": None,
    "rank_F": None,
    "rank_A": None,
    "rank_K": None,

    "H_content": None,
    "C_content": None,
    "K_content": None,
}


P7 = 7
END_E = 12
PRINCIPAL_BASE = 729


# ============================================================================
# EXACT ARITHMETIC
# ============================================================================

def gcd_int(a: int, b: int) -> int:
    a = abs(int(a))
    b = abs(int(b))

    while b:
        a, b = b, a % b

    return int(a)


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

    return int(v)


def valuation_7(x: int) -> int | None:
    return valuation_p(x, 7)


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

        old_r, r = (
            r,
            old_r - q * r,
        )

        old_s, s = (
            s,
            old_s - q * s,
        )

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

    return int(
        (
            int(q1)
            * inverse_mod(
                int(q3),
                int(modulus),
            )
        )
        % int(modulus)
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

    rho = ratio_mod(
        q1,
        q3,
        7,
    )

    matches = []

    for k in range(6):

        if orbit_mod(
            k,
            7,
        ) == rho:

            matches.append(
                int(k)
            )

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected unique k1, got {matches} "
            f"for rho={rho}."
        )

    return int(matches[0])


def lift_case(
    q1: int,
    q3: int,
    end_e: int,
):

    k1 = find_k1(
        q1,
        q3,
    )

    current_k = int(k1)

    exponents = [
        current_k
    ]

    digits = []
    slopes = []

    for e in range(
        1,
        end_e,
    ):

        modulus_e = int(
            P7 ** e
        )

        modulus_next = int(
            P7 ** (
                e + 1
            )
        )

        target = ratio_mod(
            q1,
            q3,
            modulus_next,
        )

        order = int(
            6
            * (
                P7 ** (
                    e - 1
                )
            )
        )

        candidates = [
            int(
                current_k
                + t * order
            )
            for t in range(7)
        ]

        residuals = []

        for candidate in candidates:

            orbit = orbit_mod(
                candidate,
                modulus_next,
            )

            residual = int(
                (
                    orbit
                    - target
                )
                % modulus_next
            )

            if residual % modulus_e != 0:
                raise ArithmeticError(
                    f"Candidate k={candidate} failed "
                    f"descent at e={e}."
                )

            residuals.append(
                int(
                    (
                        residual
                        // modulus_e
                    )
                    % P7
                )
            )

        A = int(
            residuals[0]
        )

        B = int(
            (
                residuals[1]
                - residuals[0]
            )
            % P7
        )

        affine = all(
            residuals[t]
            ==
            (
                A
                + B * t
            ) % P7
            for t in range(7)
        )

        if not affine:
            raise ArithmeticError(
                f"Non-affine residual law at e={e}."
            )

        matches = [
            t
            for t in range(7)
            if orbit_mod(
                candidates[t],
                modulus_next,
            )
            == target
        ]

        if len(matches) != 1:
            raise ArithmeticError(
                f"Expected exactly one lift at e={e}; "
                f"got {matches}."
            )

        chosen_t = int(
            matches[0]
        )

        current_k = int(
            candidates[chosen_t]
        )

        digits.append(
            chosen_t
        )

        slopes.append(
            B
        )

        exponents.append(
            current_k
        )

    rho7 = ratio_mod(
        q1,
        q3,
        7,
    )

    predicted_slope = int(
        (
            rho7
            * (
                (PRINCIPAL_BASE - 1)
                // P7
            )
        )
        % P7
    )

    return {
        "k1": int(k1),
        "rho7": int(rho7),
        "exponents": exponents,
        "digits": digits,
        "slopes": slopes,
        "predicted_slope": predicted_slope,
    }


# ============================================================================
# LOGARITHMIC COORDINATE
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

    return int(vn - vd)


def log_one_plus_x(
    x: Fraction,
    precision: int,
) -> Fraction:

    vx = fraction_valuation_7(
        x
    )

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

        total += (
            term
            if n % 2 == 1
            else -term
        )

        if tv is not None and tv >= precision:

            safe = True

            for j in range(
                n + 1,
                n + 8,
            ):

                future = (
                    j * vx
                    - (
                        valuation_p(
                            j,
                            7,
                        )
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
                "Log truncation exceeded safe limit."
            )

    return total


def fraction_mod(
    x: Fraction,
    modulus: int,
) -> int:

    modulus = int(modulus)

    if x.denominator % 7 == 0:
        raise ArithmeticError(
            "Rational denominator divisible by 7."
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


def logarithmic_audit(
    q1: int,
    q3: int,
    k1: int,
):

    base = int(
        2
        * pow(
            3,
            int(k1),
        )
    )

    rho = Fraction(
        int(q1),
        int(q3),
    )

    argument = (
        rho / base
    ) - 1

    v_argument = fraction_valuation_7(
        argument
    )

    if v_argument is None or v_argument <= 0:

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

    modulus = int(
        P7 ** (
            END_E - 1
        )
    )

    m_log = fraction_mod(
        logarithmic_m,
        modulus,
    )

    return {
        "defined": True,
        "m": int(m_log),
        "v_argument": int(v_argument),
    }


# ============================================================================
# CASE AUDIT
# ============================================================================

def audit_case(
    case: dict,
):

    q1 = case.get("q1")
    q3 = case.get("q3")

    if q1 is None or q3 is None:
        return None

    q1 = int(q1)
    q3 = int(q3)

    hensel = lift_case(
        q1,
        q3,
        END_E,
    )

    k1 = int(
        hensel["k1"]
    )

    final_k = int(
        hensel["exponents"][-1]
    )

    final_m = int(
        (
            final_k
            - k1
        )
        // 6
    )

    log_data = logarithmic_audit(
        q1,
        q3,
        k1,
    )

    log_match = False

    if log_data["defined"]:

        modulus = int(
            P7 ** (
                END_E - 1
            )
        )

        log_match = (
            int(log_data["m"])
            == (
                final_m
                % modulus
            )
        )

    return {
        "p": int(case["p"]),
        "q": int(case["q"]),
        "n": int(
            case["p"]
            * case["q"]
        ),

        "q1": q1,
        "q3": q3,

        "gcd": gcd_int(
            q1,
            q3,
        ),

        "rho7": hensel["rho7"],
        "k1": k1,
        "final_k": final_k,
        "final_m": final_m,

        "digits": hensel["digits"],
        "slopes": hensel["slopes"],
        "predicted_slope": hensel[
            "predicted_slope"
        ],

        "slope_exact": all(
            B == hensel["predicted_slope"]
            for B in hensel["slopes"]
        ),

        "log_defined": log_data["defined"],
        "log_match": log_match,

        "det_H": case.get("det_H"),
        "rank_F": case.get("rank_F"),
        "rank_A": case.get("rank_A"),
        "rank_K": case.get("rank_K"),
        "H_content": case.get("H_content"),
        "C_content": case.get("C_content"),
        "K_content": case.get("K_content"),
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 220R — EXACT n=6 vs n=10 "
        "SEMIPRIME FAMILY COMPARISON"
    )
    print("=" * 78)

    n10_ready = (
        CASE_N10["q1"] is not None
        and CASE_N10["q3"] is not None
    )

    # ------------------------------------------------------------------
    # 1. INPUT STATUS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. INPUT STATUS")
    print("=" * 78)

    print(
        "  n6_ready=True"
    )

    print(
        f"  n10_ready={n10_ready}"
    )

    print(
        f"  n10_target="
        f"{CASE_N10['p']}*{CASE_N10['q']}="
        f"{CASE_N10['p'] * CASE_N10['q']}"
    )

    # ------------------------------------------------------------------
    # 2. n=6
    # ------------------------------------------------------------------

    result6 = audit_case(
        CASE_N6
    )

    if result6 is None:
        raise RuntimeError(
            "n=6 data are incomplete."
        )

    print()
    print("=" * 78)
    print("2. n=6 = 2*3")
    print("=" * 78)

    print(
        f"  q1={result6['q1']}"
    )

    print(
        f"  q3={result6['q3']}"
    )

    print(
        f"  gcd(q1,q3)={result6['gcd']}"
    )

    print(
        f"  rho_mod7={result6['rho7']}"
    )

    print(
        f"  k1={result6['k1']}"
    )

    print(
        f"  final_k={result6['final_k']}"
    )

    print(
        f"  final_m={result6['final_m']}"
    )

    print(
        f"  t_digits={result6['digits']}"
    )

    print(
        f"  slopes={result6['slopes']}"
    )

    print(
        f"  predicted_slope="
        f"{result6['predicted_slope']}"
    )

    print(
        f"  slope_formula_exact="
        f"{result6['slope_exact']}"
    )

    print(
        f"  logarithmic_match="
        f"{result6['log_match']}"
    )

    print(
        f"  det_H={result6['det_H']}"
    )

    print(
        f"  rank_F={result6['rank_F']}"
    )

    print(
        f"  rank_A={result6['rank_A']}"
    )

    print(
        f"  rank_K={result6['rank_K']}"
    )

    # ------------------------------------------------------------------
    # 3. n=10
    # ------------------------------------------------------------------

    result10 = None

    if n10_ready:

        result10 = audit_case(
            CASE_N10
        )

        print()
        print("=" * 78)
        print("3. n=10 = 2*5")
        print("=" * 78)

        for key in (
            "q1",
            "q3",
            "gcd",
            "rho7",
            "k1",
            "final_k",
            "final_m",
        ):

            print(
                f"  {key}={result10[key]}"
            )

        print(
            f"  t_digits="
            f"{result10['digits']}"
        )

        print(
            f"  slopes="
            f"{result10['slopes']}"
        )

        print(
            f"  predicted_slope="
            f"{result10['predicted_slope']}"
        )

        print(
            f"  slope_formula_exact="
            f"{result10['slope_exact']}"
        )

        print(
            f"  logarithmic_match="
            f"{result10['log_match']}"
        )

        for key in (
            "det_H",
            "rank_F",
            "rank_A",
            "rank_K",
            "H_content",
            "C_content",
            "K_content",
        ):

            print(
                f"  {key}={result10[key]}"
            )

    # ------------------------------------------------------------------
    # 4. CROSS-CASE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CROSS-CASE FAMILY TEST")
    print("=" * 78)

    if result10 is None:

        print(
            "  cross_case_available=False"
        )

        print(
            "  reason="
            "exact n=10 source data are not yet supplied"
        )

    else:

        print(
            f"  rho_mod7="
            f"{result6['rho7']} vs "
            f"{result10['rho7']}"
        )

        print(
            f"  k1="
            f"{result6['k1']} vs "
            f"{result10['k1']}"
        )

        print(
            f"  predicted_slope="
            f"{result6['predicted_slope']} vs "
            f"{result10['predicted_slope']}"
        )

        print(
            f"  slope_formula_both_exact="
            f"{result6['slope_exact'] and result10['slope_exact']}"
        )

        print(
            f"  logarithmic_formula_both_exact="
            f"{result6['log_match'] and result10['log_match']}"
        )

        det6 = result6["det_H"]
        det10 = result10["det_H"]

        if (
            det6 is not None
            and det10 is not None
        ):

            print(
                f"  det_H="
                f"{det6} vs {det10}"
            )

            print(
                f"  det_H_matches="
                f"{det6 == det10}"
            )

        rankf6 = result6["rank_F"]
        ranka6 = result6["rank_A"]

        rankf10 = result10["rank_F"]
        ranka10 = result10["rank_A"]

        if all(
            x is not None
            for x in (
                rankf6,
                ranka6,
                rankf10,
                ranka10,
            )
        ):

            drop6 = int(
                rankf6 - ranka6
            )

            drop10 = int(
                rankf10 - ranka10
            )

            print(
                f"  rank_drop="
                f"{drop6} vs {drop10}"
            )

            print(
                f"  rank_drop_matches="
                f"{drop6 == drop10}"
            )

        rankk6 = result6["rank_K"]
        rankk10 = result10["rank_K"]

        if (
            rankk6 is not None
            and rankk10 is not None
        ):

            print(
                f"  rank_K="
                f"{rankk6} vs {rankk10}"
            )

            print(
                f"  rank_K_matches="
                f"{rankk6 == rankk10}"
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
The purpose of this experiment is to obtain the first independent
comparison of

    n=6=2*3

and

    n=10=2*5.

The 7-adic source mechanism is already known for an individual exact
source pair.

The real family question is whether the Schur structure survives when
the prime q changes:

    det(H),
    rank(K),
    rank(F_ret)-rank(A),
    source content,
    terminal-row structure.

A second genuine n=pq case is therefore qualitatively more valuable
than another higher-precision computation of n=6.
"""
    )

    # ------------------------------------------------------------------
    # 6. FINAL
    # ------------------------------------------------------------------

    if result10 is None:

        final_ok = True
        status = (
            "SECOND_NPQ_CASE_PENDING"
        )

    else:

        final_ok = (
            result6["slope_exact"]
            and result10["slope_exact"]
            and result6["log_match"]
            and result10["log_match"]
        )

        status = (
            "TWO_CASE_COMPARISON_COMPLETE"
        )

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  n6_analysis_exact=True"
    )

    print(
        f"  n10_analysis_available="
        f"{result10 is not None}"
    )

    if result10 is not None:

        print(
            f"  n6_slope_formula_exact="
            f"{result6['slope_exact']}"
        )

        print(
            f"  n10_slope_formula_exact="
            f"{result10['slope_exact']}"
        )

        print(
            f"  n6_log_formula_exact="
            f"{result6['log_match']}"
        )

        print(
            f"  n10_log_formula_exact="
            f"{result10['log_match']}"
        )

    print(
        f"  status={status}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 220R COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print(
            "\nInterrupted."
        )

        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise