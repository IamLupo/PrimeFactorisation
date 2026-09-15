from sympy import factorint
from math import prod


def macmahon_M1(n):
    """
    M1(n) = sigma(n), computed directly from divisors.
    """
    total = 0
    for d in range(1, n + 1):
        if n % d == 0:
            total += n // d
    return total


def squarefree_M1_from_factors(factors):
    """
    Product (p + 1), valid only when all prime exponents are 1.
    """
    return prod(p + 1 for p in factors)


def general_M1_from_factorization(factorization):
    """
    Product over p^a:
        (p^(a+1) - 1)/(p - 1)
    """
    result = 1

    for p, a in factorization.items():
        result *= (p ** (a + 1) - 1) // (p - 1)

    return result


test_numbers = [
    6,
    15,
    30,
    42,
    210,
    12,
    18,
    20,
    36,
    72,
    108,
    180,
    300,
]

print("START EXPERIMENT 1")
print()

print(
    "n\tfactorization\tM1_real\t"
    "M1_squarefree\tM1_general\t"
    "squarefree_error\tgeneral_error"
)

for n in test_numbers:
    fac = factorint(n)

    factors = list(fac.keys())

    M1_real = macmahon_M1(n)
    M1_sq = squarefree_M1_from_factors(factors)
    M1_general = general_M1_from_factorization(fac)

    print(
        f"{n}\t"
        f"{fac}\t"
        f"{M1_real}\t"
        f"{M1_sq}\t"
        f"{M1_general}\t"
        f"{M1_sq - M1_real}\t"
        f"{M1_general - M1_real}"
    )

print()
print("FINISHED EXPERIMENT 1")
