from itertools import combinations
from sympy import symbols, expand, factor, symmetrize
from sympy import primerange

from m import (
    macmahon_M1,
    macmahon_M2,
    macmahon_M3_direct,
    calc_M1,
    calc_M1_old,
    calc_M2,
    calc_x,
    get_R,
)


# =============================================================================
# CONFIGURATION
# =============================================================================

MAX_PRIME = 50
MAX_N = 2000


# =============================================================================
# HELPERS
# =============================================================================

def factor_dict(primes):
    """
    Convert a squarefree prime tuple into the factor dictionary
    expected by the new m.py functions.

    Example:
        (2, 3, 5, 7)
    becomes:
        {2: 1, 3: 1, 5: 1, 7: 1}
    """
    return {p: 1 for p in primes}


def product(values):
    result = 1

    for x in values:
        result *= x

    return result


def n_from_primes(primes):
    return product(primes)


def elementary_symmetric_4(primes):
    p, q, r, s = primes

    e1 = p + q + r + s

    e2 = (
        p*q
        + p*r
        + p*s
        + q*r
        + q*s
        + r*s
    )

    e3 = (
        p*q*r
        + p*q*s
        + p*r*s
        + q*r*s
    )

    e4 = p*q*r*s

    return e1, e2, e3, e4


# =============================================================================
# TEST 1
# M1: DIRECT VS FACTORIZED
# =============================================================================

def test_M1(prime_sets):

    print("=" * 80)
    print("TEST 1: M1 DIRECT VS FACTORIZED")
    print("=" * 80)

    failures = []

    for primes in prime_sets:

        n = n_from_primes(primes)
        factors = factor_dict(primes)

        direct = macmahon_M1(n)
        closed = calc_M1(factors)
        old = calc_M1_old(primes)

        if direct != closed or direct != old:
            failures.append(
                (
                    primes,
                    n,
                    direct,
                    closed,
                    old
                )
            )

    if not failures:
        print("PASS: all M1 formulas agree.")
    else:
        print(f"FAILURES: {len(failures)}")

        for failure in failures[:20]:
            print(failure)

    print()


# =============================================================================
# TEST 2
# M2: DIRECT VS NEW CLOSED FORM
# =============================================================================

def test_M2(prime_sets):

    print("=" * 80)
    print("TEST 2: M2 DIRECT VS NEW CLOSED FORM")
    print("=" * 80)

    failures = []

    for primes in prime_sets:

        n = n_from_primes(primes)
        factors = factor_dict(primes)

        M1 = calc_M1(factors)

        direct = macmahon_M2(n)
        closed = calc_M2(n, factors, M1)

        if direct != closed:
            failures.append(
                (
                    primes,
                    n,
                    direct,
                    closed,
                    direct - closed
                )
            )

    if not failures:
        print("PASS: all M2 formulas agree.")
    else:
        print(f"FAILURES: {len(failures)}")

        for failure in failures[:20]:
            print(failure)

    print()


# =============================================================================
# TEST 3
# R DIRECT VS STRUCTURAL PRODUCT
# =============================================================================

def test_R_structural(prime_sets):

    print("=" * 80)
    print("TEST 3: R DIRECT VS STRUCTURAL PRODUCT")
    print("=" * 80)

    failures = []

    for primes in prime_sets:

        n = n_from_primes(primes)
        factors = factor_dict(primes)

        M1 = calc_M1(factors)
        M2 = calc_M2(n, factors, M1)

        R_direct = get_R(n, M1, M2)

        x = calc_x(factors)

        R_structural = (
            M1 * (n**2 - n + 1)
            - x
        )

        if R_direct != R_structural:
            failures.append(
                (
                    primes,
                    n,
                    M1,
                    M2,
                    x,
                    R_direct,
                    R_structural,
                    R_direct - R_structural,
                )
            )

    if not failures:
        print(
            "PASS: 4-body structural identity holds "
            "for every test."
        )
    else:
        print(f"FAILURES: {len(failures)}")

        for failure in failures[:20]:
            print(failure)

    print()


# =============================================================================
# TEST 4
# X = PRODUCT(p^3 + 1)
# =============================================================================

def test_x(prime_sets):

    print("=" * 80)
    print("TEST 4: calc_x VS EXPLICIT PRODUCT")
    print("=" * 80)

    failures = []

    for primes in prime_sets:

        factors = factor_dict(primes)

        x_direct = product(
            p**3 + 1
            for p in primes
        )

        x_calc = calc_x(factors)

        if x_direct != x_calc:
            failures.append(
                (
                    primes,
                    x_direct,
                    x_calc
                )
            )

    if not failures:
        print("PASS: calc_x is exactly the cubic product.")
    else:
        print(f"FAILURES: {len(failures)}")

        for failure in failures[:20]:
            print(failure)

    print()


# =============================================================================
# TEST 5
# SYMBOLIC 4-BODY POLYNOMIAL
# =============================================================================

def derive_R4():

    print("=" * 80)
    print("TEST 5: SYMBOLIC 4-BODY R(e1,e2,e3,e4)")
    print("=" * 80)

    p, q, r, s = symbols(
        "p q r s"
    )

    variables = [
        p,
        q,
        r,
        s
    ]

    n = p*q*r*s

    M1 = (
        (p + 1)
        * (q + 1)
        * (r + 1)
        * (s + 1)
    )

    x = (
        (p**3 + 1)
        * (q**3 + 1)
        * (r**3 + 1)
        * (s**3 + 1)
    )

    R = expand(
        M1 * (n**2 - n + 1)
        - x
    )

    print("Raw R4:")
    print()
    print(R)
    print()

    result = symmetrize(
        R,
        variables,
        formal=True
    )

    symmetric_R = result[0]
    remainder = result[1]
    mapping = result[2]

    print("Symmetric R4:")
    print()
    print(factor(symmetric_R))
    print()

    print("Expanded symmetric R4:")
    print()
    print(expand(symmetric_R))
    print()

    print("Mapping:")
    print(mapping)
    print()

    print("Remainder:")
    print(remainder)
    print()

    if remainder == 0:
        print(
            "PASS: R4 is completely symmetric."
        )
    else:
        print(
            "WARNING: nonzero remainder."
        )

    print()

    return symmetric_R, mapping


# =============================================================================
# TEST 6
# VERIFY SYMBOLIC FORMULA NUMERICALLY
# =============================================================================

def test_symbolic_R4(
    symmetric_R,
    prime_sets
):

    print("=" * 80)
    print("TEST 6: SYMBOLIC R4 VS NUMERICAL R")
    print("=" * 80)

    e1, e2, e3, e4 = symbols(
        "e1 e2 e3 e4"
    )

    p, q, r, s = symbols(
        "p q r s"
    )

    R_raw = expand(
        (p + 1)
        * (q + 1)
        * (r + 1)
        * (s + 1)
        * ((p*q*r*s)**2 - p*q*r*s + 1)
        -
        (p**3 + 1)
        * (q**3 + 1)
        * (r**3 + 1)
        * (s**3 + 1)
    )

    result = symmetrize(
        R_raw,
        [p, q, r, s],
        formal=True
    )

    R_sym = result[0]

    # symmetrize normally returns s1,s2,s3,s4.
    s1, s2, s3, s4 = symbols(
        "s1 s2 s3 s4"
    )

    R_e = expand(
        R_sym.subs({
            s1: e1,
            s2: e2,
            s3: e3,
            s4: e4,
        })
    )

    print("R4(e1,e2,e3,e4):")
    print()
    print(factor(R_e))
    print()

    failures = []

    for primes in prime_sets:

        values = elementary_symmetric_4(
            primes
        )

        symbolic_value = R_e.subs({
            e1: values[0],
            e2: values[1],
            e3: values[2],
            e4: values[3],
        })

        n = n_from_primes(primes)
        factors = factor_dict(primes)

        M1 = calc_M1(factors)
        M2 = calc_M2(
            n,
            factors,
            M1
        )

        numerical_R = get_R(
            n,
            M1,
            M2
        )

        if symbolic_value != numerical_R:
            failures.append(
                (
                    primes,
                    numerical_R,
                    symbolic_value,
                    numerical_R - symbolic_value
                )
            )

    if not failures:
        print(
            "PASS: 4-body symmetric polynomial "
            "holds for every test."
        )
    else:
        print(
            f"FAILURES: {len(failures)}"
        )

        for failure in failures[:20]:
            print(failure)

    print()

    return R_e


# =============================================================================
# TEST 7
# FACTOR BY M1
# =============================================================================

def factor_R4_by_M1(R_e):

    print("=" * 80)
    print("TEST 7: DOES R4 FACTOR BY M1?")
    print("=" * 80)

    e1, e2, e3, e4 = symbols(
        "e1 e2 e3 e4"
    )

    M1_e = (
        1
        + e1
        + e2
        + e3
        + e4
    )

    print("M1(e) =")
    print(M1_e)
    print()

    quotient = factor(
        R_e / M1_e
    )

    print("R4 / M1 =")
    print()
    print(quotient)
    print()

    # Exact divisibility check.
    remainder = expand(
        R_e
    ).subs(
        e4,
        -1 - e1 - e2 - e3
    )

    if expand(remainder) == 0:
        print(
            "PASS: M1 is an exact polynomial factor."
        )
    else:
        print(
            "FAIL: M1 is not an exact factor."
        )

    print()

    return quotient


# =============================================================================
# TEST 8
# EXPLICIT EXAMPLES
# =============================================================================

def explicit_examples():

    print("=" * 80)
    print("TEST 8: EXPLICIT 4-BODY EXAMPLES")
    print("=" * 80)

    examples = [
        (2, 3, 5, 7),
        (2, 3, 5, 11),
        (2, 3, 5, 13),
        (2, 3, 7, 11),
        (2, 5, 7, 11),
        (3, 5, 7, 11),
    ]

    for primes in examples:

        n = n_from_primes(primes)
        factors = factor_dict(primes)

        M1 = calc_M1(factors)
        M2 = calc_M2(
            n,
            factors,
            M1
        )

        x = calc_x(factors)

        R = get_R(
            n,
            M1,
            M2
        )

        R_struct = (
            M1 * (n**2 - n + 1)
            - x
        )

        e1, e2, e3, e4 = (
            elementary_symmetric_4(
                primes
            )
        )

        print()
        print(f"primes = {primes}")
        print(f"n      = {n}")
        print(f"e1     = {e1}")
        print(f"e2     = {e2}")
        print(f"e3     = {e3}")
        print(f"e4     = {e4}")
        print(f"M1     = {M1}")
        print(f"M2     = {M2}")
        print(f"x      = {x}")
        print(f"R      = {R}")
        print(f"Rstruct= {R_struct}")

    print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    primes = list(
        primerange(
            2,
            MAX_PRIME + 1
        )
    )

    prime_sets = []

    for combo in combinations(
        primes,
        4
    ):

        n = n_from_primes(combo)

        if n <= MAX_N:
            prime_sets.append(combo)

    print()
    print("=" * 80)
    print("4-BODY STRUCTURAL INVESTIGATION")
    print("=" * 80)
    print()
    print(f"MAX_PRIME = {MAX_PRIME}")
    print(f"MAX_N     = {MAX_N}")
    print(f"TEST CASES = {len(prime_sets)}")
    print()

    test_M1(prime_sets)

    test_M2(prime_sets)

    test_x(prime_sets)

    test_R_structural(prime_sets)

    symmetric_R, mapping = derive_R4()

    R_e = test_symbolic_R4(
        symmetric_R,
        prime_sets
    )

    factor_R4_by_M1(
        R_e
    )

    explicit_examples()

    print("=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == "__main__":
    main()