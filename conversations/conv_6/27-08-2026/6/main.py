import sympy
import math


# ============================================================
# SETTINGS
# ============================================================

MODULI = [3, 5, 7, 11, 13, 17]


# ============================================================
# DECOMPOSE
# ============================================================

def decompose(value, modulus):
    quotient = value // modulus
    remainder = value % modulus

    return quotient, remainder


# ============================================================
# GENERATE RANDOM PRIMES
# ============================================================

p = sympy.randprime(100, 1000)
q = sympy.randprime(100, 1000)

n = p * q


# ============================================================
# HEADER
# ============================================================

print("=" * 80)
print("MULTIPLE MODULUS FACTOR DECOMPOSITION")
print("=" * 80)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()


# ============================================================
# FINGERPRINTS
# ============================================================

print("=" * 80)
print("FINGERPRINTS OF n")
print("=" * 80)

for m in MODULI:
    print(f"n mod {m:2} = {n % m:2}")

print()


# ============================================================
# DECOMPOSE p AND q FOR EVERY MODULUS
# ============================================================

print("=" * 80)
print("DECOMPOSITIONS")
print("=" * 80)

for m in MODULI:

    x, y = decompose(p, m)
    k, l = decompose(q, m)

    print()
    print(f"MODULUS m = {m}")
    print("-" * 80)

    print(f"p = {m} * {x} + {y}")
    print(f"q = {m} * {k} + {l}")

    # Expanded product
    term1 = m * m * x * k
    term2 = m * x * l
    term3 = m * k * y
    term4 = y * l

    reconstructed = term1 + term2 + term3 + term4

    print()
    print("Expansion:")
    print(
        f"{m*m}*x*k + "
        f"{m}*x*l + "
        f"{m}*k*y + "
        f"y*l"
    )

    print(
        f"= {term1} + {term2} + {term3} + {term4}"
    )

    print(f"= {reconstructed}")

    print(f"Correct: {reconstructed == n}")


# ============================================================
# COMPARE TWO MODULI
# ============================================================

print()
print("=" * 80)
print("COMPARE DIFFERENT DECOMPOSITIONS")
print("=" * 80)


for i in range(len(MODULI) - 1):

    m1 = MODULI[i]
    m2 = MODULI[i + 1]

    px, py = decompose(p, m1)
    qx, qy = decompose(q, m1)

    px2, py2 = decompose(p, m2)
    qx2, qy2 = decompose(q, m2)

    print()
    print("=" * 80)
    print(f"m1 = {m1}    m2 = {m2}")
    print("=" * 80)

    print()
    print("p:")
    print(
        f"{m1} * {px} + {py} = "
        f"{m2} * {px2} + {py2}"
    )

    print(
        f"{m1}*{px} - {m2}*{px2} = "
        f"{py2} - {py}"
    )

    left_p = m1 * px - m2 * px2
    right_p = py2 - py

    print(
        f"{left_p} = {right_p}   "
        f"Correct: {left_p == right_p}"
    )

    print()
    print("q:")
    print(
        f"{m1} * {qx} + {qy} = "
        f"{m2} * {qx2} + {qy2}"
    )

    print(
        f"{m1}*{qx} - {m2}*{qx2} = "
        f"{qy2} - {qy}"
    )

    left_q = m1 * qx - m2 * qx2
    right_q = qy2 - qy

    print(
        f"{left_q} = {right_q}   "
        f"Correct: {left_q == right_q}"
    )


# ============================================================
# PRODUCT EQUALITY BETWEEN TWO DECOMPOSITIONS
# ============================================================

print()
print("=" * 80)
print("PRODUCT REPRESENTATIONS")
print("=" * 80)

for i in range(len(MODULI) - 1):

    m1 = MODULI[i]
    m2 = MODULI[i + 1]

    x, y = decompose(p, m1)
    k, l = decompose(q, m1)

    xp, yp = decompose(p, m2)
    kp, lp = decompose(q, m2)

    left = (m1 * x + y) * (m1 * k + l)
    right = (m2 * xp + yp) * (m2 * kp + lp)

    print()
    print(f"m1 = {m1}")
    print(
        f"({m1}*{x} + {y})"
        f"({m1}*{k} + {l})"
    )

    print()

    print(f"m2 = {m2}")
    print(
        f"({m2}*{xp} + {yp})"
        f"({m2}*{kp} + {lp})"
    )

    print()
    print(f"Left  = {left}")
    print(f"Right = {right}")
    print(f"Equal = {left == right}")


# ============================================================
# FINAL FINGERPRINT
# ============================================================

print()
print("=" * 80)
print("FINAL")
print("=" * 80)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")

print()

print("Fingerprint:")
print([n % m for m in MODULI])