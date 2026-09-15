import sympy as sp
from itertools import product

# ============================================================
# SYMBOLS
# ============================================================

p, q, r, s = sp.symbols('p q r s')
e1, e2, e3, e4 = sp.symbols('e1 e2 e3 e4')

# ============================================================
# KAPPA FORMULAS
# ============================================================

k3 = (
    -e1**2
    + e1*e2
    + e1*e3
    + e1
    - e2**2
    + e2*e3
    + e2
    - e3**2
    - 2*e3
)

k4 = (
    -e1**2
    + e1*e2
    + e1*e3
    - 2*e1*e4
    + e1
    - e2**2
    + e2*e3
    + e2*e4
    + e2
    - e3**2
    + e3*e4
    - 2*e3
    - e4**2
    + e4
)

# ============================================================
# ELEMENTARY SYMMETRIC POLYNOMIALS
# ============================================================

E1_3 = p + q + r
E2_3 = p*q + p*r + q*r
E3_3 = p*q*r

E1_4 = p + q + r + s
E2_4 = p*q + p*r + p*s + q*r + q*s + r*s
E3_4 = p*q*r + p*q*s + p*r*s + q*r*s
E4_4 = p*q*r*s

K3 = sp.expand(
    k3.subs({
        e1: E1_3,
        e2: E2_3,
        e3: E3_3
    })
)

K4 = sp.expand(
    k4.subs({
        e1: E1_4,
        e2: E2_4,
        e3: E3_4,
        e4: E4_4
    })
)

# ============================================================
# BASIC TRANSITION
# ============================================================

transition = sp.factor(K4 - K3)

print("=" * 80)
print("EXACT 3 -> 4 TRANSITION")
print("=" * 80)
print()
print("K4 - K3 =")
print(transition)
print()

expected = sp.expand(s * (1 - s) * (1 - K3))

print("Verification of:")
print("K4 - K3 = s(1-s)(1-K3)")
print()
print("remainder =", sp.factor(transition - expected))
print()

# ============================================================
# MODULAR SYMBOLIC TESTS
# ============================================================

def modular_substitution_test(expr, variable, value, name):
    result = sp.factor(expr.subs(variable, value))
    print(f"{name}:")
    print(result)
    print()
    return result


print("=" * 80)
print("SPECIAL SUBSTITUTIONS FOR CONTROLLED s")
print("=" * 80)
print()

modular_substitution_test(K4 - K3, s, 0, "s = 0")
modular_substitution_test(K4 - K3, s, 1, "s = 1")
modular_substitution_test(K4 - K3, s, -1, "s = -1")

# ============================================================
# FACTORIZATION WITH RESPECT TO s
# ============================================================

print("=" * 80)
print("FACTORIZATION IN s")
print("=" * 80)
print()

poly_s = sp.Poly(K4 - K3, s)

print("Degree in s:", poly_s.degree())
print()
print("Coefficients in s:")
for power, coeff in enumerate(reversed(poly_s.all_coeffs())):
    print(f"s^{power}:")
    print(sp.factor(coeff))
    print()

print("Factor:")
print(sp.factor(K4 - K3))
print()

# ============================================================
# DIRECT DIFFERENCE FROM THE CONTROLLED-VARIABLE LAW
# ============================================================

transition_residual = sp.factor(
    (K4 - K3) - s*(1-s)*(1-K3)
)

print("=" * 80)
print("TRANSITION RESIDUAL")
print("=" * 80)
print()
print(transition_residual)
print()

# ============================================================
# EXPRESS K3 DIRECTLY IN p,q,r
# ============================================================

print("=" * 80)
print("KAPPA_3 IN p,q,r")
print("=" * 80)
print()

print(sp.factor(K3))
print()

# ============================================================
# SPECIALIZE r
# ============================================================

special_r_values = [
    0,
    1,
    -1,
    2,
    3,
    4,
    5
]

print("=" * 80)
print("KAPPA_3 UNDER CONTROLLED r")
print("=" * 80)
print()

for rv in special_r_values:
    value = sp.factor(K3.subs(r, rv))
    print(f"r = {rv}:")
    print(value)
    print()

# ============================================================
# LOOK FOR DEPENDENCE ON p+q AND pq
# ============================================================

u, v = sp.symbols('u v')

# Substitute:
# p+q = u
# pq = v
#
# Since K3 is symmetric in p,q, we can reduce it
# by symmetric reduction.

print("=" * 80)
print("KAPPA_3 AS A POLYNOMIAL IN p+q AND pq")
print("=" * 80)
print()

# Sympy symmetric reduction
sym_k3 = sp.symmetrize(
    K3,
    [p, q],
    formal=True
)

print(sym_k3)
print()

# ============================================================
# MANUAL SUBSTITUTION USING:
#
# p + q = u
# pq = v
#
# Rewrite powers using q = u-p and eliminate p.
# ============================================================

K3_u_v = sp.expand(K3.subs(q, u-p))

# Reduce modulo p^2 - u p + v
reduced = sp.rem(
    sp.Poly(K3_u_v, p),
    sp.Poly(p**2 - u*p + v, p)
)

K3_uv = sp.factor(reduced.as_expr())

print("Reduced K3:")
print(K3_uv)
print()

# ============================================================
# MODULAR COLLAPSE SEARCH
# ============================================================

print("=" * 80)
print("SYMBOLIC COLLAPSE SEARCH")
print("=" * 80)
print()

tests = {
    "K3 mod r": sp.factor(K3.subs(r, 0)),
    "K3 mod (r-1)": sp.factor(K3.subs(r, 1)),
    "K3 mod (r+1)": sp.factor(K3.subs(r, -1)),
    "K4-K3 mod s": sp.factor((K4-K3).subs(s, 0)),
    "K4-K3 mod (s-1)": sp.factor((K4-K3).subs(s, 1)),
    "K4-K3 mod (s+1)": sp.factor((K4-K3).subs(s, -1)),
}

for name, expr in tests.items():
    print(name)
    print(expr)
    print()

# ============================================================
# NUMERICAL MODULAR EXPERIMENTS
# ============================================================

def kappa_numeric(xs):
    """
    Exact definition:
        P = product(1 + x_i^3) / product(1+x_i)
        kappa = 1-P
    """
    numerator = 1
    denominator = 1

    for x in xs:
        numerator *= (1 + x**3)
        denominator *= (1 + x)

    return sp.Rational(1) - sp.Rational(numerator, denominator)


def modular_experiment(
    max_prime=30,
    max_value=10
):
    """
    Search for residue patterns involving controlled r,s.

    We compare:
        kappa3(p,q,r)
        kappa4(p,q,r,s)

    modulo small primes.
    """

    primes = list(sp.primerange(2, max_prime + 1))

    interesting = []

    for prime in primes:

        for pv in range(2, max_value + 1):
            for qv in range(2, max_value + 1):

                if pv == qv:
                    continue

                for rv in range(0, max_value + 1):

                    k3v = kappa_numeric([pv, qv, rv])

                    for sv in range(0, max_value + 1):

                        k4v = kappa_numeric([pv, qv, rv, sv])

                        # Only consider cases where denominators
                        # are invertible modulo prime.
                        den3 = (1+pv)*(1+qv)*(1+rv)
                        den4 = den3*(1+sv)

                        if den3 % prime == 0:
                            continue

                        if den4 % prime == 0:
                            continue

                        k3mod = int(k3v.p) * pow(
                            int(k3v.q), -1, prime
                        ) % prime

                        k4mod = int(k4v.p) * pow(
                            int(k4v.q), -1, prime
                        ) % prime

                        diff = (k4mod - k3mod) % prime

                        # Check controlled-variable prediction:
                        predicted = (
                            sv
                            * (1-sv)
                            * (1-k3mod)
                        ) % prime

                        if diff != predicted:
                            print(
                                "ERROR:",
                                prime,
                                pv,
                                qv,
                                rv,
                                sv
                            )
                            return

                        # Interesting if kappa loses dependence
                        # on p,q after fixing r,s.
                        interesting.append(
                            (
                                prime,
                                pv,
                                qv,
                                rv,
                                sv,
                                k3mod,
                                k4mod
                            )
                        )

    print("Numerical modular transition identity verified.")
    print("Number of valid samples:", len(interesting))

    return interesting


# ============================================================
# DEPENDENCE TEST
# ============================================================

def test_dependence_on_pq(
    prime,
    r_value,
    s_value=None,
):
    sample_max=prime

    """
    Test whether kappa3 or kappa4 becomes independent
    of p,q modulo a prime after fixing r,s.
    """

    values3 = {}

    for pv in range(1, sample_max):
        for qv in range(1, sample_max):

            if (1+pv) % prime == 0:
                continue
            if (1+qv) % prime == 0:
                continue
            if (1+r_value) % prime == 0:
                continue

            value = kappa_numeric([pv, qv, r_value])

            residue = (
                int(value.p)
                * pow(int(value.q), -1, prime)
            ) % prime

            values3[(pv, qv)] = residue

    unique3 = sorted(set(values3.values()))

    print()
    print("=" * 80)
    print("DEPENDENCE TEST")
    print("=" * 80)
    print()
    print("prime =", prime)
    print("r =", r_value)
    print("distinct kappa3 residues =", unique3)
    print("number of distinct residues =", len(unique3))

    if len(unique3) == 1:
        print(">>> POSSIBLE COLLAPSE: kappa3 independent of p,q")
    else:
        print("kappa3 still depends on p,q")

    if s_value is not None:

        values4 = {}

        for pv in range(1, sample_max):
            for qv in range(1, sample_max):

                if (1+pv) % prime == 0:
                    continue
                if (1+qv) % prime == 0:
                    continue
                if (1+r_value) % prime == 0:
                    continue
                if (1+s_value) % prime == 0:
                    continue

                value = kappa_numeric(
                    [pv, qv, r_value, s_value]
                )

                residue = (
                    int(value.p)
                    * pow(int(value.q), -1, prime)
                ) % prime

                values4[(pv, qv)] = residue

        unique4 = sorted(set(values4.values()))

        print()
        print("s =", s_value)
        print("distinct kappa4 residues =", unique4)
        print("number of distinct residues =", len(unique4))

        if len(unique4) == 1:
            print(">>> POSSIBLE COLLAPSE: kappa4 independent of p,q")
        else:
            print("kappa4 still depends on p,q")


# ============================================================
# SEARCH FOR SPECIAL r,s MODULI
# ============================================================

def search_special_controls(
    prime_limit=30,
    control_limit=10
):
    """
    Search for r,s values where kappa3 or kappa4
    exhibits unusually small residue sets.
    """

    print("=" * 80)
    print("SEARCHING FOR SPECIAL CONTROLLED VARIABLES")
    print("=" * 80)

    primes = list(sp.primerange(2, prime_limit + 1))

    discoveries = []

    for prime in primes:

        for rv in range(control_limit + 1):

            if (1 + rv) % prime == 0:
                continue

            residues3 = set()

            for pv in range(1, prime):
                for qv in range(1, prime):

                    if (1+pv) % prime == 0:
                        continue
                    if (1+qv) % prime == 0:
                        continue

                    value = kappa_numeric(
                        [pv, qv, rv]
                    )

                    residue = (
                        int(value.p)
                        * pow(int(value.q), -1, prime)
                    ) % prime

                    residues3.add(residue)

            if len(residues3) <= 2:

                discoveries.append(
                    (
                        "kappa3",
                        prime,
                        rv,
                        sorted(residues3)
                    )
                )

                print(
                    "Interesting:",
                    "prime =", prime,
                    "r =", rv,
                    "kappa3 residues =", sorted(residues3)
                )

        # Four-body controlled search
        for rv in range(control_limit + 1):

            if (1 + rv) % prime == 0:
                continue

            for sv in range(control_limit + 1):

                if (1 + sv) % prime == 0:
                    continue

                residues4 = set()

                for pv in range(1, prime):
                    for qv in range(1, prime):

                        if (1+pv) % prime == 0:
                            continue
                        if (1+qv) % prime == 0:
                            continue

                        value = kappa_numeric(
                            [pv, qv, rv, sv]
                        )

                        residue = (
                            int(value.p)
                            * pow(int(value.q), -1, prime)
                        ) % prime

                        residues4.add(residue)

                if len(residues4) <= 2:

                    discoveries.append(
                        (
                            "kappa4",
                            prime,
                            rv,
                            sv,
                            sorted(residues4)
                        )
                    )

                    print(
                        "Interesting:",
                        "prime =", prime,
                        "r =", rv,
                        "s =", sv,
                        "kappa4 residues =",
                        sorted(residues4)
                    )

    return discoveries


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    print()
    print("KAPPA CONTROLLED-VARIABLE STRUCTURE ANALYZER")
    print()

    # Exact symbolic tests
    print("Running symbolic tests...")
    print()

    # Numerical transition verification
    # Uncomment for larger experiments.
    #
    # data = modular_experiment(
    #     max_prime=29,
    #     max_value=8
    # )

    # Test selected cases
    test_dependence_on_pq(
        prime=7,
        r_value=2,
        s_value=3
    )

    test_dependence_on_pq(
        prime=11,
        r_value=2,
        s_value=3
    )

    # Search automatically for unusual modular collapse.
    #
    # discoveries = search_special_controls(
    #     prime_limit=29,
    #     control_limit=8
    # )

    print()
    print("=" * 80)
    print("DONE")
    print("=" * 80)
