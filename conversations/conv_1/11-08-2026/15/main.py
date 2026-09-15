import sympy as sp

# ============================================================
# FACTOR-RECOVERY EXPERIMENT
# ============================================================

def F(x):
    return x**2 - x + 1


def hidden_pair_invariant(p, q):
    return F(p) * F(q)


def recover_sum_candidates(n, C):
    """
    Given:
        n = pq
        C = F(p)F(q)

    Solve for u = p+q.
    """

    u = sp.symbols('u')

    equation = sp.expand(
        u**2
        - (n + 1)*u
        + (n**2 - 3*n + 1 - C)
    )

    print("Equation for p+q:")
    print(sp.factor(equation))

    roots = sp.solve(equation, u)

    print()
    print("Candidate values of p+q:")
    print(roots)

    factors = []

    for root in roots:

        if root.is_integer:
            root = int(root)

            D = root**2 - 4*n

            if D >= 0 and sp.sqrt(D).is_integer:

                d = int(sp.sqrt(D))

                if (root + d) % 2 == 0:

                    pp = (root + d)//2
                    qq = (root - d)//2

                    if pp * qq == n:
                        factors.append((pp, qq))

    print()
    print("Recovered factor pairs:")
    print(factors)

    return factors


# ============================================================
# TEST WITH KNOWN FACTORS
# ============================================================

def test_factor_recovery(p, q, r):

    n = p*q

    k3 = 1 - F(p)*F(q)*F(r)

    C = sp.Rational(
        1 - k3,
        F(r)
    )

    print("=" * 80)
    print("FACTOR RECOVERY TEST")
    print("=" * 80)

    print("p =", p)
    print("q =", q)
    print("n =", n)
    print("r =", r)

    print()
    print("kappa3 =", k3)

    print()
    print("Recovered F(p)F(q) =", C)

    print()

    return recover_sum_candidates(n, C)


# Example:
test_factor_recovery(
    p=17,
    q=23,
    r=5
)
