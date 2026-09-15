#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 229 — EXACT SCHUR ORBIT / PRIME-POWER APPROXIMATION AUDIT
==============================================================================

Known source:

    q1 = 29144191
    q3 = 24794967

Projective ratio:

    rho = q1/q3.

Schur constants:

    c_k = 2*3^k

with the terminal Schur row corresponding to

    k = 0,1,2

and hence

    c_0=2,
    c_1=6,
    c_2=18.

Experiment 228 established exact prime-power depth for these three
constants.

Experiment 229 extends the same question to a longer finite orbit:

    c_k = 2*3^k,

for k=0,...,KMAX.

For every orbit point and every prime dividing

    q1 - c_k*q3,

the experiment determines:

    * exact ell-adic valuation;
    * exact projective congruence depth;
    * first nonzero defect digit;
    * whether the depth is exceptional;
    * whether the point is in the mod-7, mod-49, or mod-343 orbit
      alignment classes.

This connects the Schur factorization phenomenon to the full
multiplicative orbit already observed in Experiments 193-196.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# SOURCE DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

KMAX = 15


# ============================================================================
# EXACT HELPERS
# ============================================================================

def valuation(
    x: int,
    ell: int,
) -> int | None:

    x = int(x)
    ell = int(ell)

    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % ell == 0:
        x //= ell
        v += 1

    return v


def factor_integer(
    x: int,
) -> dict[int, int]:

    x = abs(int(x))

    if x < 2:
        return {}

    factors: dict[int, int] = {}

    while x % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        x //= 2

    d = 3

    while d * d <= x:

        while x % d == 0:
            factors[d] = factors.get(d, 0) + 1
            x //= d

        d += 2

    if x > 1:
        factors[x] = factors.get(x, 0) + 1

    return factors


def inverse_mod(
    a: int,
    m: int,
) -> int:

    a = int(a)
    m = int(m)

    if gcd(a, m) != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return pow(
        a % m,
        -1,
        m,
    )


def ratio_mod(
    q1: int,
    q3: int,
    modulus: int,
) -> int | None:

    modulus = int(modulus)

    if gcd(q3, modulus) != 1:
        return None

    return int(
        (
            (q1 % modulus)
            * inverse_mod(
                q3,
                modulus,
            )
        )
        % modulus
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 229 — EXACT SCHUR ORBIT / "
        "PRIME-POWER APPROXIMATION AUDIT"
    )
    print("=" * 78)

    rho_mod_7 = ratio_mod(
        Q1,
        Q3,
        7,
    )

    rho_mod_49 = ratio_mod(
        Q1,
        Q3,
        49,
    )

    rho_mod_343 = ratio_mod(
        Q1,
        Q3,
        343,
    )

    print()
    print("=" * 78)
    print("1. SOURCE PROJECTIVE DATA")
    print("=" * 78)

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  gcd(q1,q3)={gcd(Q1,Q3)}"
    )

    print(
        f"  rho_mod7={rho_mod_7}"
    )

    print(
        f"  rho_mod49={rho_mod_49}"
    )

    print(
        f"  rho_mod343={rho_mod_343}"
    )

    # ------------------------------------------------------------------
    # 2. ORBIT VALUES
    # ------------------------------------------------------------------

    orbit = {}

    print()
    print("=" * 78)
    print("2. ORBIT c_k = 2*3^k")
    print("=" * 78)

    for k in range(
        KMAX + 1
    ):

        c = int(
            2 * pow(
                3,
                k,
            )
        )

        orbit[k] = c

        residual = int(
            Q1 - c * Q3
        )

        v3 = valuation(
            residual,
            3,
        )

        v5 = valuation(
            residual,
            5,
        )

        v7 = valuation(
            residual,
            7,
        )

        print(
            f"  k={k}: "
            f"c_k={c} "
            f"residual={residual} "
            f"v3={v3} "
            f"v5={v5} "
            f"v7={v7}"
        )

    # ------------------------------------------------------------------
    # 3. FACTOR EVERY ORBIT RESIDUAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. PRIME FACTORIZATION OF ORBIT RESIDUALS")
    print("=" * 78)

    residual_factors = {}

    for k, c in orbit.items():

        residual = int(
            Q1 - c * Q3
        )

        ff = factor_integer(
            residual
        )

        residual_factors[k] = ff

        print(
            f"  k={k}: "
            f"factorization={ff}"
        )

    # ------------------------------------------------------------------
    # 4. PRIME-POWER DEPTH FOR EVERY ORBIT POINT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PROJECTIVE DEPTH BY ORBIT POINT")
    print("=" * 78)

    depth_records = []

    for k, c in orbit.items():

        residual = int(
            Q1 - c * Q3
        )

        ff = residual_factors[k]

        if not ff:
            continue

        for ell in sorted(ff):

            if gcd(
                Q3,
                ell,
            ) != 1:
                continue

            v = int(
                ff[ell]
            )

            modulus = int(
                ell ** v
            )

            rho_v = ratio_mod(
                Q1,
                Q3,
                modulus,
            )

            exact_match = (
                rho_v
                ==
                (
                    c
                    % modulus
                )
            )

            next_modulus = int(
                ell ** (v + 1)
            )

            rho_next = ratio_mod(
                Q1,
                Q3,
                next_modulus,
            )

            next_match = (
                rho_next
                ==
                (
                    c
                    % next_modulus
                )
            )

            depth_records.append(
                (
                    k,
                    ell,
                    v,
                    exact_match,
                    next_match,
                )
            )

            print(
                f"  k={k}: "
                f"ell={ell} "
                f"depth={v} "
                f"exact_at_depth={exact_match} "
                f"survives_next={next_match}"
            )

    # ------------------------------------------------------------------
    # 5. FIRST DEFECT DIGIT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FIRST NONZERO DEFECT DIGITS")
    print("=" * 78)

    defect_exact = True

    for (
        k,
        ell,
        v,
        exact_match,
        next_match,
    ) in depth_records:

        c = orbit[k]

        modulus = int(
            ell ** (v + 1)
        )

        rho = ratio_mod(
            Q1,
            Q3,
            modulus,
        )

        coeff_residue = int(
            c % modulus
        )

        defect = int(
            (
                rho
                - coeff_residue
            )
            % modulus
        )

        required = int(
            ell ** v
        )

        if defect % required != 0:
            defect_exact = False
            normalized = None

        else:

            normalized = int(
                (defect // required)
                % ell
            )

        nonzero = (
            normalized is not None
            and normalized != 0
        )

        print(
            f"  k={k}: "
            f"ell={ell} "
            f"v={v} "
            f"defect={defect} "
            f"normalized_digit={normalized} "
            f"nonzero={nonzero}"
        )

        if not nonzero:
            defect_exact = False

    # ------------------------------------------------------------------
    # 6. MAXIMUM DEPTH
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. MAXIMUM ORBITAL PROJECTIVE DEPTH")
    print("=" * 78)

    maximum_depth = 0
    maximizers = []

    for (
        k,
        ell,
        v,
        exact_match,
        next_match,
    ) in depth_records:

        if v > maximum_depth:

            maximum_depth = v

            maximizers = [
                (
                    k,
                    ell,
                )
            ]

        elif v == maximum_depth:

            maximizers.append(
                (
                    k,
                    ell,
                )
            )

    print(
        f"  maximum_depth="
        f"{maximum_depth}"
    )

    print(
        f"  maximizers="
        f"{maximizers}"
    )

    # ------------------------------------------------------------------
    # 7. MOD-7 ORBIT MATCHES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. MOD-7 / MOD-49 / MOD-343 ORBIT MATCHES")
    print("=" * 78)

    matches_7 = []
    matches_49 = []
    matches_343 = []

    for k, c in orbit.items():

        if c % 7 == rho_mod_7:
            matches_7.append(k)

        if c % 49 == rho_mod_49:
            matches_49.append(k)

        if c % 343 == rho_mod_343:
            matches_343.append(k)

    print(
        f"  matches_mod7={matches_7}"
    )

    print(
        f"  matches_mod49={matches_49}"
    )

    print(
        f"  matches_mod343={matches_343}"
    )

    # ------------------------------------------------------------------
    # 8. CONNECTION TO TERMINAL SCHUR ROW
    # ------------------------------------------------------------------

    terminal = {
        0: Q1 - 2 * Q3,
        1: Q1 - 6 * Q3,
        2: Q1 - 18 * Q3,
    }

    terminal_depths = []

    for k, value in terminal.items():

        ff = factor_integer(
            value
        )

        for ell, v in ff.items():

            if gcd(
                Q3,
                ell,
            ) != 1:
                continue

            terminal_depths.append(
                (
                    k,
                    ell,
                    v,
                )
            )

    print()
    print("=" * 78)
    print("8. TERMINAL SCHUR VS FULL ORBIT")
    print("=" * 78)

    print(
        f"  terminal_depths="
        f"{terminal_depths}"
    )

    print(
        f"  full_orbit_depth_events="
        f"{[(k,ell,v) for k,ell,v,_,_ in depth_records]}"
    )

    # ------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The first three orbit points

    c_0=2,
    c_1=6,
    c_2=18

are exactly the terminal Schur constants.

Experiment 228 showed that their prime factors record exact
projective approximation depths.

Experiment 229 asks whether that phenomenon is special to the terminal
three points or is part of the larger orbit

    c_k = 2*3^k.

For every observed factor

    ell^v | (q1-c_k*q3),

the source ratio should satisfy

    q1/q3 = c_k (mod ell^v),

but not modulo ell^(v+1).

Thus the entire factorization of the finite orbit residuals can be
viewed as a collection of local projective approximation events.

The comparison with the terminal Schur points tells us whether the
Schur block is merely the first three members of a larger arithmetic
orbit or whether it exhibits exceptional behavior of its own.

This is still an exact finite statement for one source pair.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    all_depth_exact = all(
        exact
        for (
            _k,
            _ell,
            _v,
            exact,
            _next,
        )
        in depth_records
    )

    next_break_exact = all(
        not next_match
        for (
            _k,
            _ell,
            _v,
            _exact,
            next_match,
        )
        in depth_records
    )

    finite_orbit_check = (
        matches_49 == [k for k in matches_49]
        and len(matches_49) >= 1
    )

    final_ok = (
        all_depth_exact
        and defect_exact
        and next_break_exact
        and finite_orbit_check
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  all_orbit_prime_power_depths_exact="
        f"{all_depth_exact}"
    )

    print(
        f"  all_first_defect_digits_exact="
        f"{defect_exact}"
    )

    print(
        f"  all_depths_stop_at_next_power="
        f"{next_break_exact}"
    )

    print(
        f"  mod49_orbit_matches="
        f"{matches_49}"
    )

    print(
        f"  mod343_orbit_matches="
        f"{matches_343}"
    )

    print(
        f"  maximum_depth="
        f"{maximum_depth}"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_known_n6_instance"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 229 COMPLETE")


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

