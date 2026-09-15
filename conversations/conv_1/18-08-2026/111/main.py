#!/usr/bin/env python3

from fractions import Fraction
from math import factorial, gcd
from functools import reduce


# ============================================================================
# EXPERIMENT 114 — EXACT RESIDUAL-QUOTIENT FALLING-FACTORIAL AUDIT
# ============================================================================
#
# Goal:
#
#   Experiment 113 established
#
#       C_p(k) = (K-k)_s Q_p(k)
#
#   where
#
#       s = max(p-2, ceil((p-c)/2)).
#
#   Experiment 114 studies the residual polynomial Q_p(k) directly.
#
#   We expand Q_p(k) in the falling-factorial basis
#
#       Q_p(k) = sum_r q[p,r] * k_(r)
#
#   with
#
#       k_(r) = k(k-1)...(k-r+1).
#
#   Everything is exact over QQ.
#   Floating point = forbidden.
#   SymPy = not used.
#
# ============================================================================


# ----------------------------------------------------------------------------
# Exact data
# ----------------------------------------------------------------------------

A_even = {
    0: [
        Fraction(-12879, 1),
        Fraction(-28241, 1),
        Fraction(-26989, 1),
        Fraction(-13611, 1),
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
        Fraction(0, 1),
    ],
    4: [
        Fraction(-2083937, 30720),
        Fraction(-83529, 640),
        Fraction(-7965025, 129024),
        Fraction(2225141, 46080),
        Fraction(710501, 215040),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
    6: [
        Fraction(85591, 61440),
        Fraction(83651, 46080),
        Fraction(-3174439, 2580480),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
    8: [
        Fraction(-4913, 491520),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
}

A_odd = {
    1: [
        Fraction(12143, 560),
        Fraction(2699231, 1344),
        Fraction(9120441, 2240),
        Fraction(21975383, 6720),
        Fraction(4460869, 5760),
        Fraction(42929, 1680),
        Fraction(0, 1),
    ],
    3: [
        Fraction(-989, 11520),
        Fraction(-12024227, 46080),
        Fraction(-175956721, 322560),
        Fraction(-67903883, 161280),
        Fraction(-2590159, 53760),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
    5: [
        Fraction(-517, 23040),
        Fraction(396119, 30720),
        Fraction(26625517, 1290240),
        Fraction(-234707, 129024),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
    7: [
        Fraction(373, 1290240),
        Fraction(-1028053, 5160960),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
}

B_even = {
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
        Fraction(0, 1),
    ],
    4: [
        Fraction(129415, 3072),
        Fraction(267779, 3840),
        Fraction(224417, 4480),
        Fraction(-101119, 5040),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
    6: [
        Fraction(-2267, 5120),
        Fraction(-6053, 11520),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
}

B_odd = {
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
        Fraction(0, 1),
    ],
    5: [
        Fraction(4457, 46080),
        Fraction(-16819, 5760),
        Fraction(-5769, 896),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
    7: [
        Fraction(-421, 322560),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
        Fraction(0, 1),
    ],
}


CHANNELS = {
    "A-even": {"data": A_even, "K": 6, "c": 0},
    "A-odd":  {"data": A_odd,  "K": 6, "c": 0},
    "B-even": {"data": B_even, "K": 5, "c": 1},
    "B-odd":  {"data": B_odd, "K": 5, "c": 1},
}


# ----------------------------------------------------------------------------
# Basic polynomial utilities.
#
# Coefficients are stored in ascending powers:
#
#   [a0, a1, ..., an] = a0 + a1*k + ... + an*k^n
# ----------------------------------------------------------------------------

def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    if not out:
        return [Fraction(0)]
    return out


def is_zero_poly(poly):
    return all(c == 0 for c in poly)


def poly_degree(poly):
    p = trim(poly)
    if is_zero_poly(p):
        return -1
    return len(p) - 1


def poly_add(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += b[i]
    return trim(out)


def poly_sub(a, b):
    n = max(len(a), len(b))
    out = [Fraction(0)] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] -= b[i]
    return trim(out)


def poly_scale(a, c):
    return trim([x * c for x in a])


def poly_mul(a, b):
    if is_zero_poly(a) or is_zero_poly(b):
        return [Fraction(0)]
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return trim(out)


def poly_eval(poly, x):
    x = Fraction(x)
    acc = Fraction(0)
    for c in reversed(poly):
        acc = acc * x + c
    return acc


def poly_derivative(poly):
    if len(poly) <= 1:
        return [Fraction(0)]
    return trim([Fraction(i) * poly[i] for i in range(1, len(poly))])


def monomial(degree, coeff=Fraction(1)):
    if degree < 0:
        return [Fraction(0)]
    out = [Fraction(0)] * (degree + 1)
    out[degree] = Fraction(coeff)
    return out


def poly_divmod_exact(f, g):
    f = trim(f)
    g = trim(g)

    if is_zero_poly(g):
        raise ZeroDivisionError("polynomial division by zero")

    df = poly_degree(f)
    dg = poly_degree(g)

    if df < dg:
        return [Fraction(0)], f

    rem = f[:]
    q = [Fraction(0)] * (df - dg + 1)

    lead_g = g[dg]

    while poly_degree(rem) >= dg and not is_zero_poly(rem):
        dr = poly_degree(rem)
        shift = dr - dg
        coeff = rem[dr] / lead_g
        q[shift] += coeff

        sub = monomial(shift, coeff)
        sub = poly_mul(sub, g)
        rem = poly_sub(rem, sub)

    return trim(q), trim(rem)


def poly_div_exact(f, g):
    q, r = poly_divmod_exact(f, g)
    if not is_zero_poly(r):
        raise ValueError(
            f"non-exact polynomial division: remainder={poly_to_string(r)}"
        )
    return q


# ----------------------------------------------------------------------------
# Polynomial gcd over QQ, normalized monic.
# ----------------------------------------------------------------------------

def poly_monic(p):
    p = trim(p)
    if is_zero_poly(p):
        return [Fraction(0)]
    lead = p[-1]
    return trim([x / lead for x in p])


def poly_gcd(a, b):
    a = trim(a)
    b = trim(b)

    if is_zero_poly(a):
        return poly_monic(b)
    if is_zero_poly(b):
        return poly_monic(a)

    while not is_zero_poly(b):
        _, r = poly_divmod_exact(a, b)
        a, b = b, r

    return poly_monic(a)


# ----------------------------------------------------------------------------
# Falling factorial
#
# x_(r) = x(x-1)...(x-r+1)
# ----------------------------------------------------------------------------

def falling_poly(r):
    p = [Fraction(1)]
    for t in range(r):
        p = poly_mul(p, [Fraction(-t), Fraction(1)])
    return p


# ----------------------------------------------------------------------------
# Terminal factor:
#
# (K-k)_s = (K-k)(K-k-1)...(K-k-s+1)
#
# As a polynomial in k.
# ----------------------------------------------------------------------------

def terminal_factor(K, s):
    p = [Fraction(1)]
    for t in range(s):
        # K-k-t = (K-t) - k
        p = poly_mul(p, [Fraction(K - t), Fraction(-1)])
    return p


# ----------------------------------------------------------------------------
# Lagrange interpolation over QQ.
#
# Given values f(0),...,f(K), returns the unique polynomial of degree <= K.
# ----------------------------------------------------------------------------

def interpolate_values(values):
    K = len(values) - 1
    result = [Fraction(0)]

    for i, yi in enumerate(values):
        basis = [Fraction(1)]
        denom = Fraction(1)

        for j in range(K + 1):
            if j == i:
                continue

            # (k-j)
            basis = poly_mul(basis, [Fraction(-j), Fraction(1)])
            denom *= Fraction(i - j)

        basis = poly_scale(basis, yi / denom)
        result = poly_add(result, basis)

    return trim(result)


# ----------------------------------------------------------------------------
# Pretty printing.
# ----------------------------------------------------------------------------

def frac_to_str(q):
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    return f"{q.numerator}/{q.denominator}"


def poly_to_string(poly, var="k"):
    poly = trim(poly)

    if is_zero_poly(poly):
        return "0"

    pieces = []

    for power in range(poly_degree(poly), -1, -1):
        coeff = poly[power]
        if coeff == 0:
            continue

        sign = "+" if coeff > 0 else "-"
        mag = abs(coeff)

        if power == 0:
            body = frac_to_str(mag)
        elif power == 1:
            if mag == 1:
                body = var
            else:
                body = f"{frac_to_str(mag)}*{var}"
        else:
            if mag == 1:
                body = f"{var}**{power}"
            else:
                body = f"{frac_to_str(mag)}*{var}**{power}"

        if not pieces:
            if coeff < 0:
                pieces.append("-" + body)
            else:
                pieces.append(body)
        else:
            pieces.append(f" {sign} {body}")

    return "".join(pieces)


def list_to_str(values):
    return "[" + ", ".join(frac_to_str(v) for v in values) + "]"


# ----------------------------------------------------------------------------
# Primitive integer signature.
# ----------------------------------------------------------------------------

def primitive_integer_signature(poly):
    poly = trim(poly)

    if is_zero_poly(poly):
        return [0]

    denoms = [c.denominator for c in poly if c != 0]
    common_den = 1

    for d in denoms:
        common_den = common_den * d // gcd(common_den, d)

    ints = [int(c * common_den) for c in poly]

    common_gcd = reduce(gcd, [abs(x) for x in ints if x != 0], 0)

    if common_gcd == 0:
        return [0]

    ints = [x // common_gcd for x in ints]

    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-v for v in ints]
            break

    return ints


# ----------------------------------------------------------------------------
# Rank over QQ.
# ----------------------------------------------------------------------------

def matrix_rank(matrix):
    if not matrix:
        return 0

    A = [[Fraction(x) for x in row] for row in matrix]

    rows = len(A)
    cols = max(len(row) for row in A)

    for row in A:
        if len(row) < cols:
            row.extend([Fraction(0)] * (cols - len(row)))

    rank = 0
    col = 0

    while rank < rows and col < cols:
        pivot = None

        for r in range(rank, rows):
            if A[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            col += 1
            continue

        A[rank], A[pivot] = A[pivot], A[rank]

        pivot_value = A[rank][col]
        A[rank] = [x / pivot_value for x in A[rank]]

        for r in range(rows):
            if r == rank:
                continue
            factor = A[r][col]
            if factor != 0:
                A[r] = [
                    A[r][c] - factor * A[rank][c]
                    for c in range(cols)
                ]

        rank += 1
        col += 1

    return rank


# ----------------------------------------------------------------------------
# Support law from Experiment 111.
# ----------------------------------------------------------------------------

def support_s(p, c):
    half = (p - c + 1) // 2  # floor((p+1-c)/2)
    linear = p - 2
    return max(linear, half)


# ----------------------------------------------------------------------------
# Extract Q_p(k):
#
#   C_p(k) = (K-k)_s Q_p(k)
# ----------------------------------------------------------------------------

def extract_terminal_quotient(values, K, s):
    C = interpolate_values(values)
    divisor = terminal_factor(K, s)
    Q = poly_div_exact(C, divisor)

    # Pointwise reconstruction.
    for k in range(K + 1):
        lhs = poly_eval(C, k)
        rhs = poly_eval(divisor, k) * poly_eval(Q, k)
        if lhs != rhs:
            raise AssertionError(
                f"pointwise reconstruction failed at k={k}"
            )

    return C, divisor, Q


# ----------------------------------------------------------------------------
# Falling-factorial basis transform.
#
# If
#
#   Q(k) = sum_r q_r k_(r),
#
# recover q_r by triangular evaluation:
#
#   Q(n) = sum_{r=0}^n q_r n_(r).
# ----------------------------------------------------------------------------

def falling_basis_coefficients(poly):
    d = poly_degree(poly)

    if d < 0:
        return []

    q = [Fraction(0)] * (d + 1)

    for n in range(d + 1):
        value = poly_eval(poly, n)

        used = Fraction(0)
        for r in range(n):
            used += q[r] * factorial(n) // factorial(n - r)

        q[n] = value - used
        q[n] /= Fraction(factorial(n))

    return trim(q)


# ----------------------------------------------------------------------------
# Reconstruct from falling basis.
# ----------------------------------------------------------------------------

def falling_basis_reconstruct(coeffs):
    result = [Fraction(0)]

    for r, coeff in enumerate(coeffs):
        if coeff == 0:
            continue
        result = poly_add(result, poly_scale(falling_poly(r), coeff))

    return trim(result)


# ----------------------------------------------------------------------------
# Tests on the falling-basis quotient rows.
# ----------------------------------------------------------------------------

def first_nonzero_index(values):
    for i, v in enumerate(values):
        if v != 0:
            return i
    return None


def last_nonzero_index(values):
    for i in range(len(values) - 1, -1, -1):
        if values[i] != 0:
            return i
    return None


def coefficient_matrix(rows):
    width = max((len(r) for r in rows), default=0)
    out = []
    for row in rows:
        out.append(row + [Fraction(0)] * (width - len(row)))
    return out


def normalized_row(row):
    nz = next((x for x in row if x != 0), None)
    if nz is None:
        return row[:]
    return [x / nz for x in row]


def rank1_matrix(matrix):
    if not matrix:
        return False

    base = None

    for row in matrix:
        if any(x != 0 for x in row):
            base = row
            break

    if base is None:
        return True

    pivot_index = next(i for i, x in enumerate(base) if x != 0)

    for row in matrix:
        pivot = row[pivot_index]
        if pivot == 0:
            if any(x != 0 for x in row):
                return False
            continue

        scale = pivot / base[pivot_index]

        for j in range(len(base)):
            if row[j] != scale * base[j]:
                return False

    return True


# ----------------------------------------------------------------------------
# Search for a second falling-factorial factor in Q(k).
#
# We test k_(r), (k-a)_(r), and (K-k)_r with small r.
# ----------------------------------------------------------------------------

def candidate_factor_search(Q):
    d = poly_degree(Q)
    if d < 1:
        return []

    hits = []

    # Lower-end factors k_(r).
    for r in range(1, d + 1):
        divisor = falling_poly(r)
        q, rem = poly_divmod_exact(Q, divisor)
        if is_zero_poly(rem):
            hits.append(("k_falling", r, q))

    # Shifted lower-end factors (k-a)_(r).
    for a in range(-3, 4):
        base = [Fraction(-a), Fraction(1)]
        for r in range(1, d + 1):
            divisor = [Fraction(1)]
            for t in range(r):
                divisor = poly_mul(
                    divisor,
                    [Fraction(-(a + t)), Fraction(1)]
                )

            q, rem = poly_divmod_exact(Q, divisor)
            if is_zero_poly(rem):
                hits.append(("shifted_k", a, r, q))

    return hits


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 114 — EXACT RESIDUAL-QUOTIENT")
    print("FALLING-FACTORIAL / SECOND-LAYER AUDIT")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------------
    # 1. Support law validation
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. EXACT SUPPORT LAW")
    print("=" * 78)

    support_ok = True

    for name, cfg in CHANNELS.items():
        print(f"{name}: K={cfg['K']} c={cfg['c']}")

        for p, values in cfg["data"].items():
            observed = 0
            while (
                observed < cfg["K"]
                and values[cfg["K"] - observed] == 0
            ):
                observed += 1

            predicted = support_s(p, cfg["c"])
            exact = observed == predicted

            if not exact:
                support_ok = False

            print(
                f"  p={p}: observed_s={observed} "
                f"predicted_s={predicted} exact={exact}"
            )

        print()

    # ------------------------------------------------------------------------
    # 2. Extract residual quotient polynomials
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("2. EXACT TERMINAL-FACTOR / QUOTIENT EXTRACTION")
    print("=" * 78)

    all_quotients = {}
    all_falling_rows = {}

    quotient_ok = True

    for name, cfg in CHANNELS.items():
        print(name)

        K = cfg["K"]
        c = cfg["c"]

        all_quotients[name] = {}
        all_falling_rows[name] = {}

        for p, values in cfg["data"].items():
            s = support_s(p, c)

            try:
                C, divisor, Q = extract_terminal_quotient(
                    values, K, s
                )
                exact = True
            except Exception as exc:
                exact = False
                quotient_ok = False
                C = interpolate_values(values)
                divisor = terminal_factor(K, s)
                Q = [Fraction(0)]
                print(f"  p={p}: extraction FAILED: {exc}")

            if exact:
                qdeg = poly_degree(Q)

                print(
                    f"  p={p}: s={s} "
                    f"deg(C)={poly_degree(C)} "
                    f"deg(divisor)={poly_degree(divisor)} "
                    f"deg(Q)={qdeg} exact={exact}"
                )
                print(f"    divisor = {poly_to_string(divisor)}")
                print(f"    Q(k)    = {poly_to_string(Q)}")

            all_quotients[name][p] = Q

        print()

    # ------------------------------------------------------------------------
    # 3. Residual degree law
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("3. RESIDUAL DEGREE LAW")
    print("=" * 78)

    residual_degree_ok = True

    for name, cfg in CHANNELS.items():
        K = cfg["K"]
        c = cfg["c"]

        print(name)

        for p in cfg["data"]:
            s = support_s(p, c)
            Q = all_quotients[name][p]

            observed = poly_degree(Q)
            predicted = K - s
            exact = observed == predicted

            if not exact:
                residual_degree_ok = False

            print(
                f"  p={p}: s={s} observed={observed} "
                f"predicted={predicted} exact={exact}"
            )

        print()

    # ------------------------------------------------------------------------
    # 4. Falling-factorial basis of Q
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("4. RESIDUAL QUOTIENT FALLING-FACTORIAL BASIS")
    print("=" * 78)

    basis_ok = True

    for name, cfg in CHANNELS.items():
        print(name)

        rows = []

        for p in cfg["data"]:
            Q = all_quotients[name][p]
            coeffs = falling_basis_coefficients(Q)
            recon = falling_basis_reconstruct(coeffs)

            exact = recon == trim(Q)

            if not exact:
                basis_ok = False

            rows.append(coeffs)

            print(f"  p={p}:")
            print(f"    q[p,r] = {list_to_str(coeffs)}")
            print(f"    reconstruction = {exact}")

        all_falling_rows[name] = rows
        print()

    # ------------------------------------------------------------------------
    # 5. Falling-basis support / triangularity
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("5. SECOND-LAYER FALLING-BASIS SUPPORT")
    print("=" * 78)

    triangular_ok = True

    for name, cfg in CHANNELS.items():
        print(name)

        rows = all_falling_rows[name]

        for p, coeffs in zip(cfg["data"], rows):
            first = first_nonzero_index(coeffs)
            last = last_nonzero_index(coeffs)

            print(
                f"  p={p}: first_nonzero_r={first} "
                f"last_nonzero_r={last} "
                f"degree={poly_degree(all_quotients[name][p])}"
            )

        print()

    # ------------------------------------------------------------------------
    # 6. Primitive integer signatures
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("6. PRIMITIVE INTEGER SIGNATURES OF Q FALLING-BASIS ROWS")
    print("=" * 78)

    for name, rows in all_falling_rows.items():
        print(name)
        for p, row in zip(CHANNELS[name]["data"], rows):
            print(
                f"  p={p}: {primitive_integer_signature(row)}"
            )
        print()

    # ------------------------------------------------------------------------
    # 7. Quotient coefficient-matrix rank
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("7. SECOND-LAYER COEFFICIENT-MATRIX RANK")
    print("=" * 78)

    rank_data = {}

    for name, rows in all_falling_rows.items():
        M = coefficient_matrix(rows)
        rank = matrix_rank(M)
        r1 = rank1_matrix(M)

        rank_data[name] = rank

        print(
            f"{name}: shape=({len(M)}, {len(M[0]) if M else 0}) "
            f"rank={rank} rank-1={r1}"
        )

    print()

    # ------------------------------------------------------------------------
    # 8. Terminal-zero structure in the SECOND basis
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("8. SECOND-LAYER TERMINAL FACTOR SEARCH")
    print("=" * 78)

    second_factor_ok = True

    for name, rows in all_falling_rows.items():
        print(name)

        for p, row in zip(CHANNELS[name]["data"], rows):
            Q = all_quotients[name][p]
            hits = candidate_factor_search(Q)

            # This is diagnostic only: no assumption that hits must exist.
            if hits:
                descriptions = []

                for hit in hits:
                    if hit[0] == "k_falling":
                        descriptions.append(
                            f"k_({hit[1]})"
                        )
                    else:
                        descriptions.append(
                            f"(k-{hit[1]})_({hit[2]})"
                        )

                print(
                    f"  p={p}: exact candidate factors={descriptions}"
                )
            else:
                print(f"  p={p}: exact candidate factors=[]")

        print()

    # ------------------------------------------------------------------------
    # 9. Cross-parity comparison of quotient falling rows
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("9. CROSS-PARITY QUOTIENT FALLING-BASIS COMPARISON")
    print("=" * 78)

    paired = [
        (0, "A-even", "B-even"),
        (2, "A-even", "B-even"),
        (4, "A-even", "B-even"),
        (6, "A-even", "B-even"),
        (1, "A-odd", "B-odd"),
        (3, "A-odd", "B-odd"),
        (5, "A-odd", "B-odd"),
        (7, "A-odd", "B-odd"),
    ]

    for p, left, right in paired:
        if p not in all_quotients[left] or p not in all_quotients[right]:
            continue

        A = all_falling_rows[left][list(CHANNELS[left]["data"]).index(p)]
        B = all_falling_rows[right][list(CHANNELS[right]["data"]).index(p)]

        width = max(len(A), len(B))
        A2 = A + [Fraction(0)] * (width - len(A))
        B2 = B + [Fraction(0)] * (width - len(B))

        ratios = []
        for r in range(width):
            if B2[r] != 0:
                ratios.append((r, A2[r] / B2[r]))
            elif A2[r] != 0:
                ratios.append((r, None))

        print(f"p={p}: {left} vs {right}")
        print(f"  entrywise ratios={ratios}")
        print()

    # ------------------------------------------------------------------------
    # 10. Exact reconstruction from the complete second-layer basis
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("10. EXACT SECOND-LAYER RECONSTRUCTION")
    print("=" * 78)

    reconstruction_ok = True

    for name, rows in all_falling_rows.items():
        print(name)

        for p, row in zip(CHANNELS[name]["data"], rows):
            reconstructed = falling_basis_reconstruct(row)
            original = all_quotients[name][p]

            exact = reconstructed == trim(original)

            if not exact:
                reconstruction_ok = False

            print(f"  p={p}: exact={exact}")

        print()

    # ------------------------------------------------------------------------
    # 11. Pointwise sanity check
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("11. POINTWISE SANITY CHECK")
    print("=" * 78)

    pointwise_ok = True

    for name, cfg in CHANNELS.items():
        K = cfg["K"]
        c = cfg["c"]

        print(name)

        for p, values in cfg["data"].items():
            s = support_s(p, c)
            divisor = terminal_factor(K, s)
            Q = all_quotients[name][p]

            exact = True

            for k in range(K + 1):
                lhs = values[k]
                rhs = poly_eval(divisor, k) * poly_eval(Q, k)

                if lhs != rhs:
                    exact = False
                    break

            if not exact:
                pointwise_ok = False

            print(f"  p={p}: exact={exact}")

        print()

    # ------------------------------------------------------------------------
    # 12. Structural summary
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("12. STRUCTURAL INTERPRETATION")
    print("=" * 78)
    print()
    print("Experiment 113 established the exact first residual factor:")
    print()
    print("    C_p(k) = (K-k)_s Q_p(k)")
    print()
    print("with")
    print()
    print("    s = max(p-2, ceil((p-c)/2)).")
    print()
    print("Experiment 114 now changes basis on Q_p(k):")
    print()
    print("    Q_p(k) = sum_r q[p,r] k_(r).")
    print()
    print("The purpose is to determine whether the residual kernel")
    print("contains a SECOND factorial/combinatorial layer.")
    print()
    print("A useful positive signal would be:")
    print()
    print("    * exact support in r;")
    print("    * a new triangular pattern;");
    print("    * simple terminal or lower-end factors;");
    print("    * low matrix rank;");
    print("    * a simple A/B relation.")
    print()
    print("A negative result would mean that the two endpoint factors")
    print("already expose the main simple structure available in the")
    print("observed finite data.")
    print()
    print("Everything is exact over QQ.")
    print("No floating point.")
    print("No SymPy.")
    print("No extrapolation.")
    print("No recurrence search.")
    print()

    # ------------------------------------------------------------------------
    # 13. Final exactness
    # ------------------------------------------------------------------------

    failures = 0

    checks = {
        "support_law": support_ok,
        "terminal_quotient": quotient_ok,
        "residual_degree": residual_degree_ok,
        "falling_basis": basis_ok,
        "reconstruction": reconstruction_ok,
        "pointwise": pointwise_ok,
    }

    for ok in checks.values():
        if not ok:
            failures += 1

    print("=" * 78)
    print("13. FINAL EXACTNESS")
    print("=" * 78)

    print(f"  support_law = {support_ok}")
    print(f"  terminal_quotient = {quotient_ok}")
    print(f"  residual_degree = {residual_degree_ok}")
    print(f"  second_falling_basis = {basis_ok}")
    print(f"  second_layer_reconstruction = {reconstruction_ok}")
    print(f"  pointwise_exactness = {pointwise_ok}")
    print(f"  failures = {failures}")
    print(
        "  ALL BASIC CHECKS PASS = "
        f"{failures == 0}"
    )
    print()
    print("EXPERIMENT 114 COMPLETE")


if __name__ == "__main__":
    main()

