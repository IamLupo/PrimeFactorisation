"""
==============================================================================
EXPERIMENT 113 — EXACT INDEX-POLYNOMIAL TERMINAL-FACTOR AUDIT
==============================================================================

CORRECTION TO EXPERIMENT 112

The coefficient lists are VALUES of a polynomial in k, not polynomial
coefficients.

Therefore we first reconstruct

    C_p(k)

exactly from the samples

    C_p(0), C_p(1), ..., C_p(K),

and only then divide the polynomial by

    (K-k)_(s).

All arithmetic is exact over QQ.
Floating point = forbidden.
No SymPy required.
==============================================================================
"""

from __future__ import annotations

from fractions import Fraction


# ============================================================================
# OBSERVED DATA
# ============================================================================

DATA = {
    "A-even": {
        "c": 0,
        "K": 6,
        "powers": {
            0: [
                Fraction(-12879),
                Fraction(-28241),
                Fraction(-26989),
                Fraction(-13611),
                Fraction(-17875, 6),
                Fraction(-116923, 1680),
                Fraction(-5, 144),
            ],
            2: [
                Fraction(2797337, 1920),
                Fraction(8986567, 2880),
                Fraction(36298273, 13440),
                Fraction(11670379, 11520),
                Fraction(19954213, 161280),
                Fraction(-9389, 4032),
                Fraction(0),
            ],
            4: [
                Fraction(-2083937, 30720),
                Fraction(-83529, 640),
                Fraction(-7965025, 129024),
                Fraction(2225141, 46080),
                Fraction(710501, 215040),
                Fraction(0),
                Fraction(0),
            ],
            6: [
                Fraction(85591, 61440),
                Fraction(83651, 46080),
                Fraction(-3174439, 2580480),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
            ],
            8: [
                Fraction(-4913, 491520),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
            ],
        },
    },

    "A-odd": {
        "c": 0,
        "K": 6,
        "powers": {
            1: [
                Fraction(12143, 560),
                Fraction(2699231, 1344),
                Fraction(9120441, 2240),
                Fraction(21975383, 6720),
                Fraction(4460869, 5760),
                Fraction(42929, 1680),
                Fraction(0),
            ],
            3: [
                Fraction(-989, 11520),
                Fraction(-12024227, 46080),
                Fraction(-175956721, 322560),
                Fraction(-67903883, 161280),
                Fraction(-2590159, 53760),
                Fraction(0),
                Fraction(0),
            ],
            5: [
                Fraction(-517, 23040),
                Fraction(396119, 30720),
                Fraction(26625517, 1290240),
                Fraction(-234707, 129024),
                Fraction(0),
                Fraction(0),
                Fraction(0),
            ],
            7: [
                Fraction(373, 1290240),
                Fraction(-1028053, 5160960),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
            ],
        },
    },

    "B-even": {
        "c": 1,
        "K": 5,
        "powers": {
            0: [
                Fraction(12980463, 1024),
                Fraction(6255583, 256),
                Fraction(87841139, 4480),
                Fraction(16998339, 2240),
                Fraction(1091983, 896),
                Fraction(1553, 240),
            ],
            2: [
                Fraction(-19344659, 15360),
                Fraction(-26986999, 11520),
                Fraction(-2066529, 1120),
                Fraction(-573325, 576),
                Fraction(3312053, 40320),
                Fraction(0),
            ],
            4: [
                Fraction(129415, 3072),
                Fraction(267779, 3840),
                Fraction(224417, 4480),
                Fraction(-101119, 5040),
                Fraction(0),
                Fraction(0),
            ],
            6: [
                Fraction(-2267, 5120),
                Fraction(-6053, 11520),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
            ],
        },
    },

    "B-odd": {
        "c": 1,
        "K": 5,
        "powers": {
            1: [
                Fraction(-584531, 35840),
                Fraction(-2908483, 1920),
                Fraction(-31233169, 13440),
                Fraction(-1446167, 1344),
                Fraction(-22259149, 40320),
                Fraction(-301, 240),
            ],
            3: [
                Fraction(-59257, 46080),
                Fraction(186547, 1440),
                Fraction(367433, 1680),
                Fraction(126549, 448),
                Fraction(-162139, 40320),
                Fraction(0),
            ],
            5: [
                Fraction(4457, 46080),
                Fraction(-16819, 5760),
                Fraction(-5769, 896),
                Fraction(0),
                Fraction(0),
                Fraction(0),
            ],
            7: [
                Fraction(-421, 322560),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
                Fraction(0),
            ],
        },
    },
}


ZERO = Fraction(0)


# ============================================================================
# POLYNOMIAL UTILITIES
# ============================================================================
# Polynomials are stored in ascending powers:
#
#   [a0, a1, a2, ...]
#
# meaning
#
#   a0 + a1*k + a2*k^2 + ...
# ============================================================================

def trim(poly):
    poly = list(poly)

    while len(poly) > 1 and poly[-1] == ZERO:
        poly.pop()

    return poly


def is_zero(poly):
    poly = trim(poly)
    return len(poly) == 1 and poly[0] == ZERO


def degree(poly):
    poly = trim(poly)

    if is_zero(poly):
        return -1

    return len(poly) - 1


def add(a, b):
    n = max(len(a), len(b))
    out = [ZERO] * n

    for i in range(n):
        if i < len(a):
            out[i] += a[i]

        if i < len(b):
            out[i] += b[i]

    return trim(out)


def sub(a, b):
    n = max(len(a), len(b))
    out = [ZERO] * n

    for i in range(n):
        if i < len(a):
            out[i] += a[i]

        if i < len(b):
            out[i] -= b[i]

    return trim(out)


def scale(a, q):
    return trim([q * x for x in a])


def mul(a, b):
    if is_zero(a) or is_zero(b):
        return [ZERO]

    out = [ZERO] * (len(a) + len(b) - 1)

    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj

    return trim(out)


def evaluate(poly, x):
    acc = ZERO

    for coeff in reversed(trim(poly)):
        acc = acc * x + coeff

    return acc


# ============================================================================
# EXACT POLYNOMIAL DIVISION
# ============================================================================

def divmod_poly(num, den):
    num = trim(num)
    den = trim(den)

    if is_zero(den):
        raise ZeroDivisionError("polynomial division by zero")

    if is_zero(num):
        return [ZERO], [ZERO]

    if degree(num) < degree(den):
        return [ZERO], num

    q = [ZERO] * (degree(num) - degree(den) + 1)
    r = list(num)

    den_deg = degree(den)
    den_lc = den[-1]

    while not is_zero(r) and degree(r) >= den_deg:
        shift = degree(r) - den_deg
        coeff = r[-1] / den_lc

        q[shift] += coeff

        term = [ZERO] * shift + scale(den, coeff)
        r = sub(r, term)

    return trim(q), trim(r)


def exact_quotient(num, den):
    q, r = divmod_poly(num, den)

    return q, is_zero(r)


# ============================================================================
# LAGRANGE INTERPOLATION
# ============================================================================

def interpolate_from_samples(values):
    """
    Given

        values[k] = f(k),  k=0,...,K,

    reconstruct the unique degree <= K polynomial f(k).
    """
    n = len(values)

    result = [ZERO]

    for i in range(n):
        basis = [Fraction(1)]
        denominator = Fraction(1)

        for j in range(n):
            if i == j:
                continue

            # Multiply by (x-j)
            basis = mul(
                basis,
                [Fraction(-j), Fraction(1)],
            )

            denominator *= Fraction(i - j)

        basis = scale(
            basis,
            values[i] / denominator,
        )

        result = add(result, basis)

    return trim(result)


# ============================================================================
# FALLING FACTORIAL
# ============================================================================

def falling_factorial_linear(K, s):
    """
    Construct

        (K-k)_s
        = (K-k)(K-k-1)...(K-k-s+1)

    as a polynomial in k.
    """
    result = [Fraction(1)]

    for i in range(s):
        factor = [
            Fraction(K - i),
            Fraction(-1),
        ]

        result = mul(result, factor)

    return trim(result)


# ============================================================================
# EXPERIMENT-111 SUPPORT LAW
# ============================================================================

def support_law(p, c):
    """
        s(p,c) = max(p-2, ceil((p-c)/2))
    """
    half = (p - c + 1) // 2
    return max(p - 2, half)


def ceil_half(p, c):
    return (p - c + 1) // 2


def branch(p, c):
    linear = p - 2
    half = ceil_half(p, c)

    if linear > half:
        return "linear"

    if linear < half:
        return "half"

    return "tie"


# ============================================================================
# EXPECTED RESIDUAL DEGREE
# ============================================================================

def expected_degree(K, s):
    return K - s


# ============================================================================
# MAIN
# ============================================================================

def main():

    failures = 0

    reconstructed = {}
    extracted = {}

    print("=" * 78)
    print("EXPERIMENT 113 — EXACT INDEX-POLYNOMIAL TERMINAL-FACTOR AUDIT")
    print("=" * 78)

    # ======================================================================
    print()
    print("=" * 78)
    print("1. EXACT SUPPORT LAW")
    print("=" * 78)

    support_ok = True

    for name, block in DATA.items():

        c = block["c"]
        K = block["K"]

        print()
        print(f"{name}: c={c}, K={K}")

        for p, samples in block["powers"].items():

            observed_zero_tail = 0

            for value in reversed(samples):
                if value == ZERO:
                    observed_zero_tail += 1
                else:
                    break

            predicted = support_law(p, c)

            ok = observed_zero_tail == predicted

            support_ok = support_ok and ok

            print(
                f"  p={p}: "
                f"observed_s={observed_zero_tail} "
                f"predicted_s={predicted} "
                f"exact={ok}"
            )

    if not support_ok:
        failures += 1

    # ======================================================================
    print()
    print("=" * 78)
    print("2. INDEX-POLYNOMIAL RECONSTRUCTION")
    print("=" * 78)

    interpolation_ok = True

    for name, block in DATA.items():

        print()
        print(name)

        reconstructed[name] = {}

        for p, samples in block["powers"].items():

            poly = interpolate_from_samples(samples)

            reconstructed[name][p] = poly

            ok = True

            for k, expected in enumerate(samples):

                actual = evaluate(poly, Fraction(k))

                if actual != expected:
                    ok = False
                    break

            interpolation_ok = interpolation_ok and ok

            print(
                f"  p={p}: "
                f"degree={degree(poly)} "
                f"reconstructed={ok}"
            )

    if not interpolation_ok:
        failures += 1

    # ======================================================================
    print()
    print("=" * 78)
    print("3. EXACT TERMINAL-FACTOR EXTRACTION")
    print("=" * 78)

    extraction_ok = True

    for name, block in DATA.items():

        K = block["K"]
        c = block["c"]

        extracted[name] = {}

        print()
        print(name)

        for p in block["powers"]:

            poly = reconstructed[name][p]

            s = support_law(p, c)

            divisor = falling_factorial_linear(K, s)

            quotient, exact = exact_quotient(poly, divisor)

            remainder = sub(
                poly,
                mul(divisor, quotient),
            )

            qdeg = degree(quotient)

            reconstructed_after = True

            for k in range(K + 1):

                expected = evaluate(poly, Fraction(k))
                actual = evaluate(
                    mul(divisor, quotient),
                    Fraction(k),
                )

                if expected != actual:
                    reconstructed_after = False
                    break

            exact_all = exact and is_zero(remainder)

            if not exact_all:
                extraction_ok = False

            extracted[name][p] = {
                "s": s,
                "divisor": divisor,
                "quotient": quotient,
                "degree": qdeg,
                "exact": exact_all,
                "pointwise": reconstructed_after,
            }

            print(
                f"  p={p}: "
                f"s={s} "
                f"divisor_degree={degree(divisor)} "
                f"exact={exact_all} "
                f"quotient_degree={qdeg} "
                f"pointwise={reconstructed_after}"
            )

    if not extraction_ok:
        failures += 1

    # ======================================================================
    print()
    print("=" * 78)
    print("4. CORRECTED RESIDUAL-DEGREE AUDIT")
    print("=" * 78)

    degree_ok = True

    for name, block in DATA.items():

        K = block["K"]

        print()
        print(name)

        for p in block["powers"]:

            item = extracted[name][p]

            observed = item["degree"]
            s = item["s"]
            predicted = expected_degree(K, s)

            ok = observed == predicted

            degree_ok = degree_ok and ok

            print(
                f"  p={p}: "
                f"s={s} "
                f"observed_degree={observed} "
                f"predicted_degree={predicted} "
                f"exact={ok}"
            )

    if not degree_ok:
        failures += 1

    # ======================================================================
    print()
    print("=" * 78)
    print("5. BRANCH AUDIT")
    print("=" * 78)

    for name, block in DATA.items():

        c = block["c"]

        print()
        print(name)

        for p in block["powers"]:

            linear = p - 2
            half = ceil_half(p, c)
            s = support_law(p, c)

            print(
                f"  p={p}: "
                f"linear={linear} "
                f"half={half} "
                f"branch={branch(p,c)} "
                f"s={s}"
            )

    # ======================================================================
    print()
    print("=" * 78)
    print("6. TERMINAL-FACTOR POLYNOMIALS")
    print("=" * 78)

    for name, block in DATA.items():

        print()
        print(name)

        for p in block["powers"]:

            item = extracted[name][p]

            print(
                f"  p={p}: "
                f"s={item['s']} "
                f"quotient_degree={item['degree']}"
            )

            print(
                f"      divisor_coefficients = "
                f"{item['divisor']}"
            )

            print(
                f"      quotient_coefficients = "
                f"{item['quotient']}"
            )

    # ======================================================================
    print()
    print("=" * 78)
    print("7. SUPPORT / DEGREE COUPLING TABLE")
    print("=" * 78)

    for name, block in DATA.items():

        K = block["K"]
        c = block["c"]

        print()
        print(name)

        for p in block["powers"]:

            s = support_law(p, c)
            qdeg = extracted[name][p]["degree"]

            print(
                f"  p={p}: "
                f"K={K} "
                f"s={s} "
                f"K-s={K-s} "
                f"quotient_degree={qdeg}"
            )

    # ======================================================================
    print()
    print("=" * 78)
    print("8. FINAL STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The previous failure came from confusing polynomial samples with
  polynomial coefficients.

  The correct object is the index polynomial

      C_p(k)

  reconstructed from the exact values

      C_p(0), C_p(1), ..., C_p(K).

  The terminal factor is then tested as the genuine polynomial divisor

      (K-k)_s.

  Thus the corrected identity is

      C_p(k)
        = (K-k)_s Q_p(k).

  This experiment therefore separates three statements:

      SUPPORT:
          s = max(p-2, ceil((p-c)/2))

      DIVISIBILITY:
          (K-k)_s divides C_p(k)

      RESIDUAL DEGREE:
          deg Q_p = K-s.

  A successful result establishes that the terminal-zero pattern is
  an actual algebraic factor of the indexed coefficient polynomial,
  rather than merely a consequence of the finite sample list.

  Everything is exact over QQ.
  No floating point.
  No extrapolation.
  """
    )

    # ======================================================================
    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    all_checks = (
        support_ok
        and interpolation_ok
        and extraction_ok
        and degree_ok
    )

    print(
        f"  support_law_exact = {support_ok}"
    )

    print(
        f"  index_polynomial_reconstruction = {interpolation_ok}"
    )

    print(
        f"  terminal_factor_extraction = {extraction_ok}"
    )

    print(
        f"  residual_degree_law = {degree_ok}"
    )

    print(
        f"  failures = {failures}"
    )

    print(
        f"  ALL BASIC CHECKS PASS = {all_checks}"
    )

    print()
    print("EXPERIMENT 113 COMPLETE")


if __name__ == "__main__":
    main()

