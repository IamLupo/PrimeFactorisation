from sympy import factorint, primerange
from math import isqrt


def quotient_state(p, q, r1, r2):
    """
    Exact quotient/remainder state:
        p = a + k*r1
        q = b + l*r2
    """

    a = p % r1
    k = p // r1

    b = q % r2
    l = q // r2

    n = p * q

    T = n // (r1 * r2)

    K = k * l
    E = T - K

    return {
        "n": n,
        "p": p,
        "q": q,
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "T": T,
        "K": K,
        "E": E,
    }


def nearby_semiprime_state(n, x, r1, r2):
    """
    Factor n+x completely and inspect every factor pair.
    We only use this for experimental data generation.
    """

    nx = n + x

    if nx <= 1:
        return []

    fac = factorint(nx)

    result = []

    # Generate divisor pairs p_x * q_x = nx.
    divisors = []

    def generate_divisors(items, index=0, current=1):
        if index == len(items):
            divisors.append(current)
            return

        p, exponent = items[index]

        value = 1

        for _ in range(exponent + 1):
            generate_divisors(
                items,
                index + 1,
                current * value
            )
            value *= p

    generate_divisors(list(fac.items()))

    for px in sorted(divisors):
        qx = nx // px

        # Avoid duplicate unordered pair.
        if px > qx:
            continue

        state = quotient_state(px, qx, r1, r2)

        state["x"] = x
        state["factorization_nx"] = fac
        state["px"] = px
        state["qx"] = qx

        result.append(state)

    return result


print("START EXPERIMENT 4")
print()

# Distinct odd-prime semiprimes.
prime_list = list(primerange(3, 100))

# Fixed moduli for the experiment.
r1 = 11
r2 = 13

# Example semiprimes.
examples = [
    (17, 43),
    (19, 47),
    (23, 53),
    (29, 59),
    (31, 67),
]

# Nearby displacements.
xs = list(range(-25, 26))

for p, q in examples:

    n = p * q

    original = quotient_state(p, q, r1, r2)

    print("--------------------------------------------------")
    print(f"ORIGINAL n={n}, p={p}, q={q}")
    print(
        f"original: "
        f"k={original['k']} "
        f"l={original['l']} "
        f"K={original['K']} "
        f"T={original['T']} "
        f"E={original['E']}"
    )
    print()

    print(
        "x\tn+x\tpx\tqx\t"
        "kx\tlx\tKx\tTx\tEx\t"
        "same_T\tKx-K"
    )

    for x in xs:

        states = nearby_semiprime_state(
            n,
            x,
            r1,
            r2
        )

        for state in states:

            same_T = state["T"] == original["T"]

            print(
                f"{x}\t"
                f"{n+x}\t"
                f"{state['px']}\t"
                f"{state['qx']}\t"
                f"{state['k']}\t"
                f"{state['l']}\t"
                f"{state['K']}\t"
                f"{state['T']}\t"
                f"{state['E']}\t"
                f"{same_T}\t"
                f"{state['K'] - original['K']}"
            )

print()
print("FINISHED EXPERIMENT 4")
