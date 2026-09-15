from sympy import primerange, Poly, symbols


def elementary_symmetric(values):
    """
    Return e_0 ... e_k by iterative multiplication:
        prod(1 + x_i t)
    """
    coeffs = [1]

    for x in values:
        new = [0] * (len(coeffs) + 1)

        for i, c in enumerate(coeffs):
            new[i] += c
            new[i + 1] += c * x

        coeffs = new

    return coeffs


print("START EXPERIMENT 3")
print()

primes = list(primerange(2, 50))

# -------------------------
# 3-body
# -------------------------

print("=== 3-BODY ===")
print(
    "p\tq\tr\tn\t"
    "e1\te2\te3\t"
    "M1\tM1_from_e\t"
    "poly"
)

count = 0

for i, p in enumerate(primes):
    for j in range(i + 1, len(primes)):
        q = primes[j]

        for k in range(j + 1, len(primes)):
            r = primes[k]

            n = p * q * r

            if n > 5000:
                continue

            e = elementary_symmetric([p, q, r])

            e0, e1, e2, e3 = e

            M1 = (p + 1) * (q + 1) * (r + 1)
            M1_from_e = e0 + e1 + e2 + e3

            poly = f"x^3 - ({e1})x^2 + ({e2})x - ({e3})"

            print(
                f"{p}\t{q}\t{r}\t{n}\t"
                f"{e1}\t{e2}\t{e3}\t"
                f"{M1}\t{M1_from_e}\t"
                f"{poly}"
            )

            count += 1

            if count >= 20:
                break

        if count >= 20:
            break

    if count >= 20:
        break


# -------------------------
# 4-body
# -------------------------

print()
print("=== 4-BODY ===")
print(
    "p\tq\tr\ts\tn\t"
    "e1\te2\te3\te4\t"
    "M1\tM1_from_e\t"
    "poly"
)

count = 0

for i, p in enumerate(primes):
    for j in range(i + 1, len(primes)):
        q = primes[j]

        for k in range(j + 1, len(primes)):
            r = primes[k]

            for l in range(k + 1, len(primes)):
                s = primes[l]

                n = p * q * r * s

                if n > 10000:
                    continue

                e = elementary_symmetric([p, q, r, s])

                e0, e1, e2, e3, e4 = e

                M1 = (
                    (p + 1)
                    * (q + 1)
                    * (r + 1)
                    * (s + 1)
                )

                M1_from_e = sum(e)

                poly = (
                    f"x^4 - ({e1})x^3 + "
                    f"({e2})x^2 - ({e3})x + ({e4})"
                )

                print(
                    f"{p}\t{q}\t{r}\t{s}\t{n}\t"
                    f"{e1}\t{e2}\t{e3}\t{e4}\t"
                    f"{M1}\t{M1_from_e}\t"
                    f"{poly}"
                )

                count += 1

                if count >= 20:
                    break

            if count >= 20:
                break

        if count >= 20:
            break

    if count >= 20:
        break

print()
print("FINISHED EXPERIMENT 3")
