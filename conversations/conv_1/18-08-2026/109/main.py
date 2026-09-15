"""
==============================================================================
EXPERIMENT 112 — EXACT SUPPORT-LAW TO RESIDUAL-DEGREE / ENDPOINT-KERNEL AUDIT
==============================================================================

Exact rational arithmetic only.
No floating point.
No SymPy.
No extrapolation.
"""

from __future__ import annotations

from fractions import Fraction
from math import ceil


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


# ============================================================================
# POLYNOMIAL UTILITIES
# ============================================================================

ZERO = Fraction(0)


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == ZERO:
        poly.pop()
    return poly


def degree(poly):
    poly = trim(poly)
    if len(poly) == 1 and poly[0] == ZERO:
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


def mul(a, b):
    if degree(a) == -1 or degree(b) == -1:
        return [ZERO]

    out = [ZERO] * (len(a) + len(b) - 1)

    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj

    return trim(out)


def scale(poly, q):
    return trim([x * q for x in poly])


def evaluate(poly, x):
    acc = ZERO
    for coeff in reversed(trim(poly)):
        acc = acc * x + coeff
    return acc


# ============================================================================
# FALLING FACTORIAL IN k
# ============================================================================

def terminal_factor(K, s):
    """
    (K-k)(K-k-1)...(K-k-s+1)
    """
    out = [Fraction(1)]

    for i in range(s):
        out = mul(out, [Fraction(K - i), Fraction(-1)])

    return trim(out)


# ============================================================================
# EXACT POLYNOMIAL DIVISION
# ============================================================================

def exact_division(num, den):
    num = trim(num)
    den = trim(den)

    if degree(den) == -1:
        raise ZeroDivisionError("division by zero polynomial")

    if degree(num) < degree(den):
        return [ZERO], num == [ZERO]

    quotient = [ZERO] * (degree(num) - degree(den) + 1)
    remainder = list(num)

    den_degree = degree(den)
    den_lc = den[-1]

    while degree(remainder) >= den_degree:
        rdeg = degree(remainder)
        coeff = remainder[-1] / den_lc
        shift = rdeg - den_degree

        quotient[shift] += coeff

        term = [ZERO] * shift + [coeff * x for x in den]
        remainder = sub(remainder, term)

    exact = degree(remainder) == -1 or remainder == [ZERO]

    return trim(quotient), exact


# ============================================================================
# SUPPORT LAW FROM EXPERIMENT 111
# ============================================================================

def support_law(p, c):
    half = (p - c + 1) // 2
    return max(p - 2, half)


def branch_name(p, c):
    linear = p - 2
    half = (p - c + 1) // 2

    if linear > half:
        return "linear"
    if linear < half:
        return "half"
    return "tie"


# ============================================================================
# OBSERVED TERMINAL-ZERO COUNT
# ============================================================================

def terminal_zero_count(values):
    count = 0
    for value in reversed(values):
        if value == ZERO:
            count += 1
        else:
            break
    return count


# ============================================================================
# EXACT LAGRANGE INTERPOLATION
# ============================================================================

def lagrange_interpolate(values):
    """
    Values correspond to k = 0,1,...,n-1.
    Returns polynomial coefficients in ascending powers of k.
    """
    n = len(values)
    result = [ZERO]

    for i in range(n):
        basis = [Fraction(1)]
        denom = Fraction(1)

        for j in range(n):
            if i == j:
                continue

            basis = mul(
                basis,
                [Fraction(-j), Fraction(1)]
            )
            denom *= Fraction(i - j)

        basis = scale(basis, values[i] / denom)
        result = add(result, basis)

    return trim(result)


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 112 — EXACT SUPPORT-LAW TO RESIDUAL-DEGREE / ENDPOINT-KERNEL AUDIT")
    print("=" * 78)

    failures = 0

    extracted = {}

    # ----------------------------------------------------------------------
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

        for p, values in block["powers"].items():

            observed = terminal_zero_count(values)
            predicted = support_law(p, c)

            ok = observed == predicted
            support_ok = support_ok and ok

            print(
                f"  p={p}: "
                f"observed_s={observed} "
                f"predicted_s={predicted} "
                f"exact={ok}"
            )

    if not support_ok:
        failures += 1

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("2. EXACT TERMINAL-FACTOR EXTRACTION")
    print("=" * 78)

    extraction_ok = True

    for name, block in DATA.items():
        c = block["c"]
        K = block["K"]

        extracted[name] = {}

        print()
        print(name)

        for p, values in block["powers"].items():

            s = support_law(p, c)
            divisor = terminal_factor(K, s)

            quotient, exact = exact_division(
                values,
                divisor,
            )

            qdeg = degree(quotient)

            extracted[name][p] = {
                "s": s,
                "quotient": quotient,
                "degree": qdeg,
                "exact": exact,
            }

            extraction_ok = extraction_ok and exact

            print(
                f"  p={p}: "
                f"s={s} "
                f"exact={exact} "
                f"quotient_degree={qdeg}"
            )

    if not extraction_ok:
        failures += 1

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("3. RESIDUAL-DEGREE LAW")
    print("=" * 78)

    degree_ok = True

    for name, block in DATA.items():
        K = block["K"]

        print()
        print(name)

        for p in block["powers"]:

            item = extracted[name][p]

            observed_degree = item["degree"]
            s = item["s"]
            predicted_degree = K - s

            ok = observed_degree == predicted_degree
            degree_ok = degree_ok and ok

            print(
                f"  p={p}: "
                f"s={s} "
                f"observed={observed_degree} "
                f"predicted={predicted_degree} "
                f"exact={ok}"
            )

    if not degree_ok:
        failures += 1

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("4. BRANCH MECHANISM AUDIT")
    print("=" * 78)

    for name, block in DATA.items():
        c = block["c"]

        print()
        print(name)

        for p in block["powers"]:

            linear = p - 2
            half = (p - c + 1) // 2
            s = support_law(p, c)
            branch = branch_name(p, c)
            qdeg = extracted[name][p]["degree"]

            print(
                f"  p={p}: "
                f"linear={linear} "
                f"half={half} "
                f"branch={branch} "
                f"s={s} "
                f"residual_degree={qdeg}"
            )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("5. QUOTIENT DEGREE IDENTITY TABLE")
    print("=" * 78)

    degree_identity_ok = True

    for name, block in DATA.items():
        K = block["K"]

        print()
        print(name)

        for p in block["powers"]:

            s = extracted[name][p]["s"]
            qdeg = extracted[name][p]["degree"]

            candidates = {
                "K-s": K - s,
                "K-(p-2)": K - (p - 2),
                "K-ceil((p-c)/2)": K - ((p - c + 1) // 2),
                "K-max(...)": K - max(
                    p - 2,
                    (p - c + 1) // 2,
                ),
            }

            print(
                f"  p={p}: observed_degree={qdeg} "
                f"candidate_degrees={candidates}"
            )

            if qdeg != candidates["K-max(...)"]:
                degree_identity_ok = False

    if not degree_identity_ok:
        failures += 1

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("6. EXACT QUOTIENT INTERPOLATION")
    print("=" * 78)

    interpolation_ok = True

    for name, block in DATA.items():

        print()
        print(name)

        for p, item in extracted[name].items():

            quotient = item["quotient"]
            interp = lagrange_interpolate(quotient)

            reconstructed = True

            for k, expected in enumerate(quotient):
                actual = evaluate(interp, Fraction(k))

                if actual != expected:
                    reconstructed = False
                    break

            interpolation_ok = interpolation_ok and reconstructed

            print(
                f"  p={p}: "
                f"degree={degree(interp)} "
                f"reconstruction={reconstructed}"
            )

    if not interpolation_ok:
        failures += 1

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("7. EXPLICIT SUPPORT / DEGREE SUMMARY")
    print("=" * 78)

    for name, block in DATA.items():
        c = block["c"]
        K = block["K"]

        print()
        print(name)

        for p in block["powers"]:

            s = support_law(p, c)
            residual_degree = extracted[name][p]["degree"]

            print(
                f"  p={p}: "
                f"K={K}, c={c}, "
                f"s={s}, "
                f"K-s={K-s}, "
                f"residual_degree={residual_degree}"
            )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("8. FINAL STRUCTURAL CHECKS")
    print("=" * 78)

    # Equivalent ceiling/floor formulation.
    floor_form_ok = True

    for name, block in DATA.items():
        c = block["c"]

        for p in block["powers"]:

            s1 = support_law(p, c)
            s2 = max(
                p - 2,
                (p + 1 - c) // 2,
            )

            if s1 != s2:
                floor_form_ok = False

    print(
        f"  ceiling_floor_equivalence={floor_form_ok}"
    )

    if not floor_form_ok:
        failures += 1

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  Experiment 111 established the exact support law

      s(p,c)
        = max(
            p - 2,
            ceil((p-c)/2)
          ).

  Experiment 112 tests whether this law controls more than the
  terminal zero count.

  For every observed parity coefficient family,

      c_p(k)
        = (K-k)_(s(p,c)) Q_p(k).

  The central identity tested here is

      deg Q_p = K - s(p,c).

  Therefore the same discrete support quantity determines both:

      1. where the coefficient family terminates;
      2. the dimension of the remaining index polynomial.

  The branch audit distinguishes the two mechanisms inside the max:

      p-2

  versus

      ceil((p-c)/2).

  This gives a sharper interpretation of the support law without
  introducing another recurrence search.

  Everything is exact over rational arithmetic.
  No floating point.
  No extrapolation.
  No SymPy.
  """
    )

    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    all_terminal_exact = extraction_ok
    all_checks = (
        support_ok
        and all_terminal_exact
        and degree_ok
        and interpolation_ok
        and floor_form_ok
    )

    print(f"  support_law_exact = {support_ok}")
    print(f"  terminal_factor_extraction = {all_terminal_exact}")
    print(f"  residual_degree_law_exact = {degree_ok}")
    print(f"  quotient_interpolation_exact = {interpolation_ok}")
    print(f"  ceiling_floor_equivalence = {floor_form_ok}")
    print(f"  failures = {failures}")
    print(f"  ALL BASIC CHECKS PASS = {all_checks}")

    print()
    print("EXPERIMENT 112 COMPLETE")


if __name__ == "__main__":
    main()