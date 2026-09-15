from sympy import primerange, symbols, simplify, factor
from math import isqrt


def macmahon_M1(n):
    total = 0

    for s in range(1, n + 1):
        if n % s == 0:
            total += n // s

    return total


def macmahon_M2(n):
    total = 0

    for s1 in range(1, n):
        for s2 in range(s1 + 1, n + 1):

            for m1 in range(1, n // s1 + 1):
                remainder = n - m1 * s1

                if remainder > 0 and remainder % s2 == 0:
                    m2 = remainder // s2
                    total += m1 * m2

    return total


def R_detector(n, M1, M2):
    return M1 * (n**2 - 3*n + 2) - 8*M2


print("START EXPERIMENT 2")
print()

primes = list(primerange(2, 100))

print(
    "n\tp\tq\tS\t"
    "R_real\tR_formula\t"
    "P2\tP2_formula\t"
    "S_recovered\t"
    "p_recovered\tq_recovered\tOK"
)

for i, p in enumerate(primes):
    for q in primes[i + 1:]:

        n = p * q

        if n > 1000:
            continue

        S = p + q

        M1 = macmahon_M1(n)
        M2 = macmahon_M2(n)

        R = R_detector(n, M1, M2)

        R_formula = (p**2 - 1) * (q**2 - 1) * (p + q)

        P2 = (1 - p + p**2) * (1 - q + q**2)

        P2_formula = (
            S**2
            - (n + 1) * S
            + (n - 1)**2
        )

        # From
        #
        # P2 = S^2 - (n+1)S + (n-1)^2
        #
        # solve quadratic:
        #
        # S^2 - (n+1)S + ((n-1)^2 - P2) = 0
        #
        discriminant = (
            (n + 1)**2
            - 4 * ((n - 1)**2 - P2)
        )

        sqrt_disc = isqrt(discriminant)

        recovered_candidates = []

        if sqrt_disc * sqrt_disc == discriminant:
            recovered_candidates = [
                ((n + 1) + sqrt_disc) // 2,
                ((n + 1) - sqrt_disc) // 2,
            ]

        S_recovered = None
        p_recovered = None
        q_recovered = None

        for candidate_S in recovered_candidates:
            if candidate_S == S:

                factor_disc = candidate_S**2 - 4*n

                if factor_disc >= 0:
                    d = isqrt(factor_disc)

                    if d*d == factor_disc:
                        p_recovered = (candidate_S - d) // 2
                        q_recovered = (candidate_S + d) // 2
                        S_recovered = candidate_S

        ok = (
            R == R_formula
            and P2 == P2_formula
            and S_recovered == S
            and set([p_recovered, q_recovered]) == set([p, q])
        )

        print(
            f"{n}\t{p}\t{q}\t{S}\t"
            f"{R}\t{R_formula}\t"
            f"{P2}\t{P2_formula}\t"
            f"{S_recovered}\t"
            f"{p_recovered}\t{q_recovered}\t"
            f"{ok}"
        )

print()
print("FINISHED EXPERIMENT 2")
