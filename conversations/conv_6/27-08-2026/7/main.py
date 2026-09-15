import sympy
import math


# ============================================================
# MODULUS SETS
# ============================================================

MODULI1 = [3, 5, 7, 11, 13]
MODULI2 = [17, 19, 23, 29, 31]

ALL_MODULI = MODULI1 + MODULI2


# ============================================================
# RANDOM PRIME FACTORS
# ============================================================

p = sympy.randprime(1000, 10000)
q = sympy.randprime(1000, 10000)

n = p * q


# ============================================================
# DECOMPOSITION
# ============================================================

def decompose(value, m):
    x = value // m
    r = value % m
    return x, r


# ============================================================
# PRODUCT TERMS
# ============================================================

def product_terms(p, q, m):
    x, r = decompose(p, m)
    y, s = decompose(q, m)

    A = m * m * x * y
    B = m * (x * s + y * r)
    C = r * s

    total = A + B + C

    return {
        "x": x,
        "r": r,
        "y": y,
        "s": s,
        "A": A,
        "B": B,
        "C": C,
        "total": total
    }


# ============================================================
# HEADER
# ============================================================

print("=" * 80)
print("MULTI-FINGERPRINT QUOTIENT/REMAINDER EXPERIMENT")
print("=" * 80)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()


# ============================================================
# FINGERPRINTS
# ============================================================

fp1 = [n % m for m in MODULI1]
fp2 = [n % m for m in MODULI2]

print("=" * 80)
print("FINGERPRINT 1")
print("=" * 80)

print(f"MODULI1 = {MODULI1}")
print(f"F1(n)   = {fp1}")
print()

print("=" * 80)
print("FINGERPRINT 2")
print("=" * 80)

print(f"MODULI2 = {MODULI2}")
print(f"F2(n)   = {fp2}")
print()


# ============================================================
# DECOMPOSITIONS
# ============================================================

records = {}

print("=" * 80)
print("DECOMPOSITIONS")
print("=" * 80)

for m in ALL_MODULI:

    d = product_terms(p, q, m)
    records[m] = d

    print()
    print(f"m = {m}")
    print("-" * 80)

    print(f"p = {m} * {d['x']} + {d['r']}")
    print(f"q = {m} * {d['y']} + {d['s']}")

    print()

    print(f"A = m²xy                 = {d['A']}")
    print(f"B = m(xs + yr)           = {d['B']}")
    print(f"C = rs                   = {d['C']}")
    print(f"A + B + C                = {d['total']}")

    print(f"Correct                  = {d['total'] == n}")


# ============================================================
# TEST BASIC INVARIANTS
# ============================================================

print()
print("=" * 80)
print("INVARIANT TESTS")
print("=" * 80)

for m in ALL_MODULI:

    d = records[m]

    x = d["x"]
    r = d["r"]
    y = d["y"]
    s = d["s"]

    # n mod m
    residue_product = (r * s) % m

    # quotient after removing the residue contribution
    quotient_part = (n - r * s) // m

    expected_quotient_part = m * x * y + x * s + y * r

    print()
    print(f"m = {m}")

    print(f"n mod m                 = {n % m}")
    print(f"(r*s) mod m             = {residue_product}")

    print(f"(n-r*s)/m               = {quotient_part}")
    print(
        f"m*x*y + x*s + y*r      = "
        f"{expected_quotient_part}"
    )

    print(
        f"Quotient identity       = "
        f"{quotient_part == expected_quotient_part}"
    )


# ============================================================
# COMPARE MODULI
# ============================================================

print()
print("=" * 80)
print("CROSS-MODULUS COMPARISON")
print("=" * 80)

for i in range(len(ALL_MODULI)):

    m1 = ALL_MODULI[i]

    for j in range(i + 1, len(ALL_MODULI)):

        m2 = ALL_MODULI[j]

        a = records[m1]
        b = records[m2]

        # ----------------------------------------------------
        # p relationship
        # ----------------------------------------------------

        p_left = m1 * a["x"] - m2 * b["x"]
        p_right = b["r"] - a["r"]

        # ----------------------------------------------------
        # q relationship
        # ----------------------------------------------------

        q_left = m1 * a["y"] - m2 * b["y"]
        q_right = b["s"] - a["s"]

        # ----------------------------------------------------
        # Combined quantity
        # ----------------------------------------------------

        residue_difference = (
            b["r"] * b["s"]
            - a["r"] * a["s"]
        )

        print()
        print(f"{m1} -> {m2}")

        print(
            f"p: {m1}*{a['x']} - "
            f"{m2}*{b['x']} = "
            f"{p_left}"
        )

        print(
            f"   r2-r1 = {p_right}"
        )

        print(
            f"   valid = {p_left == p_right}"
        )

        print(
            f"q: {m1}*{a['y']} - "
            f"{m2}*{b['y']} = "
            f"{q_left}"
        )

        print(
            f"   s2-s1 = {q_right}"
        )

        print(
            f"   valid = {q_left == q_right}"
        )

        print(
            f"C2-C1 = "
            f"(r2*s2)-(r1*s1) = "
            f"{residue_difference}"
        )


# ============================================================
# NEW PART:
# SEARCH FOR A SIMPLE CROSS-MODULUS RELATION
# ============================================================

print()
print("=" * 80)
print("CROSS-MODULUS NORMALIZED VALUES")
print("=" * 80)

for i in range(len(ALL_MODULI)):

    m1 = ALL_MODULI[i]

    for j in range(i + 1, len(ALL_MODULI)):

        m2 = ALL_MODULI[j]

        a = records[m1]
        b = records[m2]

        # Difference in quotient terms
        dx = m1 * a["x"] - m2 * b["x"]
        dy = m1 * a["y"] - m2 * b["y"]

        # These must equal remainder differences
        dr = b["r"] - a["r"]
        ds = b["s"] - a["s"]

        # A combined expression
        cross1 = dx * dy
        cross2 = dr * ds

        print()
        print(f"m1={m1:2}, m2={m2:2}")

        print(f"dx = {dx}")
        print(f"dr = {dr}")

        print(f"dy = {dy}")
        print(f"ds = {ds}")

        print(f"dx*dy = {cross1}")
        print(f"dr*ds = {cross2}")


# ============================================================
# COMPARE THE TWO FINGERPRINT GROUPS
# ============================================================

print()
print("=" * 80)
print("GROUP COMPARISON")
print("=" * 80)

group1_values = []
group2_values = []

for m in MODULI1:
    d = records[m]

    group1_values.append(
        (
            m,
            d["r"],
            d["s"],
            d["A"],
            d["B"],
            d["C"]
        )
    )

for m in MODULI2:
    d = records[m]

    group2_values.append(
        (
            m,
            d["r"],
            d["s"],
            d["A"],
            d["B"],
            d["C"]
        )
    )


print()
print("GROUP 1:")

for item in group1_values:
    print(item)


print()
print("GROUP 2:")

for item in group2_values:
    print(item)


# ============================================================
# FINAL
# ============================================================

print()
print("=" * 80)
print("FINAL")
print("=" * 80)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")

print()
print(f"F1 = {fp1}")
print(f"F2 = {fp2}")
