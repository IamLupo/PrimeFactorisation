#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 513
==============================================================================
CENTERED SPECTRAL POLYNOMIAL / HISTORICAL CONIC FACTORIZATION
==============================================================================
"""

import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

U, V, r = sp.symbols("U V r", nonzero=True)
N, S, Delta = sp.symbols("N S Delta")
z, xi = sp.symbols("z xi")


# ============================================================================
# HELPERS
# ============================================================================

def ex(expr):
    return sp.factor(sp.expand(expr))


def zero(expr):
    return sp.expand(expr) == 0


def certificate(label, expr):
    e = ex(expr)
    print(f"  {label}")
    print(f"    difference = {e}")
    print(f"    PASS = {e == 0}")
    print()


# ============================================================================
# START
# ============================================================================

print("=" * 78)
print("EXPERIMENT 513 START")
print("=" * 78)
print("CENTERED SPECTRAL POLYNOMIAL / HISTORICAL CONIC FACTORIZATION")
print("=" * 78)


# ============================================================================
# [1] HISTORICAL HYPERBOLIC COORDINATES
# ============================================================================

print()
print("[1] HISTORICAL HYPERBOLIC COORDINATES")
print("-" * 78)

A0 = U - V
B0 = U + V

K = ex(U**2 - V**2)

print("  A0 = U-V")
print("  B0 = U+V")
print()
print(f"  A0+B0 = {ex(A0 + B0)}")
print(f"  B0-A0 = {ex(B0 - A0)}")
print(f"  A0*B0 = {ex(A0 * B0)}")
print()

certificate(
    "A0*B0 = U^2-V^2",
    A0 * B0 - (U**2 - V**2),
)


# ============================================================================
# [2] NORMALIZED KAPPA BASES
# ============================================================================

print()
print("[2] FOUR KAPPA BASES")
print("-" * 78)

P = ex((U - V) / r)
Q = ex((U + V) / r)

bases = [
    P,
    Q,
    P + 1,
    Q + 1,
]

for i, base in enumerate(bases, start=1):
    print(f"  lambda_{i} = {base}")

print()


# ============================================================================
# [3] FOUR-BASE CHARACTERISTIC POLYNOMIAL
# ============================================================================

print()
print("[3] RAW CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi = ex(sp.prod(z - b for b in bases))

print("  chi(z) =")
print(f"    {chi}")
print()


# ============================================================================
# [4] FACTORIZATION INTO TWO QUADRATICS
# ============================================================================

print()
print("[4] TWO-QUADRATIC FACTORIZATION")
print("-" * 78)

quadratic_0 = ex(
    (z - P) * (z - Q)
)

quadratic_1 = ex(
    (z - P - 1) * (z - Q - 1)
)

print("  Q0(z) =")
print(f"    {quadratic_0}")

print()
print("  Q1(z) =")
print(f"    {quadratic_1}")

print()

certificate(
    "chi(z)=Q0(z)*Q1(z)",
    chi - quadratic_0 * quadratic_1,
)


# ============================================================================
# [5] EXPRESS THE TWO QUADRATICS THROUGH N,S
# ============================================================================

print()
print("[5] MODERN SYMMETRIC FORM")
print("-" * 78)

N_modern = ex(P * Q)
S_modern = ex(P + Q)

q0_target = z**2 - S_modern * z + N_modern
q1_target = (
    z**2
    - (S_modern + 2) * z
    + (N_modern + S_modern + 1)
)

certificate(
    "Q0 = z^2-Sz+N",
    quadratic_0 - q0_target,
)

certificate(
    "Q1 = z^2-(S+2)z+(N+S+1)",
    quadratic_1 - q1_target,
)


# ============================================================================
# [6] CENTER OF THE FOUR ROOTS
# ============================================================================

print()
print("[6] SPECTRAL CENTER")
print("-" * 78)

root_sum = ex(sum(bases))
center = ex(root_sum / 4)

print(f"  sum(lambda_i) = {root_sum}")
print(f"  spectral center = {center}")
print()

center_expected = ex(S_modern / 2 + sp.Rational(1, 2))

certificate(
    "center=(S+1)/2",
    center - center_expected,
)


# ============================================================================
# [7] CENTERED ROOT COORDINATES
# ============================================================================

print()
print("[7] CENTERED ROOTS")
print("-" * 78)

centered_roots = [
    ex(base - center)
    for base in bases
]

for i, root in enumerate(centered_roots, start=1):
    print(f"  centered lambda_{i} = {root}")

print()


# ============================================================================
# [8] CENTERED CHARACTERISTIC POLYNOMIAL
# ============================================================================

print()
print("[8] CENTERED CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

Delta_modern = ex((P - Q)**2)

chi_centered = ex(
    chi.subs(z, xi + center)
)

print("  chi(xi + (S+1)/2) =")
print(f"    {chi_centered}")
print()


# ============================================================================
# [9] CLOSED CENTERED FORM
# ============================================================================

print()
print("[9] CLOSED CENTERED FORM")
print("-" * 78)

centered_target = ex(
    (
        xi**2
        + sp.Rational(1, 4)
        - Delta_modern / 4
    )**2
    - xi**2
)

print("  proposed:")
print("    (xi^2 + 1/4 - Delta/4)^2 - xi^2")
print()

certificate(
    "centered characteristic polynomial",
    chi_centered - centered_target,
)


# ============================================================================
# [10] CENTERED FACTORIZATION
# ============================================================================

print()
print("[10] CENTERED FACTORIZATION")
print("-" * 78)

factorized_centered = ex(
    (
        (xi + sp.Rational(1, 2))**2
        - Delta_modern / 4
    )
    *
    (
        (xi - sp.Rational(1, 2))**2
        - Delta_modern / 4
    )
)

print("  factorized form:")
print(f"    {factorized_centered}")
print()

certificate(
    "centered factorization",
    chi_centered - factorized_centered,
)


# ============================================================================
# [11] HISTORICAL U,V SPECTRAL FORM
# ============================================================================

print()
print("[11] HISTORICAL U,V SPECTRAL FORM")
print("-" * 78)

spectral_poly_w = ex(
    (z - (U - V))
    * (z - (U + V))
    * (z - (U - V + r))
    * (z - (U + V + r))
)

print("  R(z) =")
print(f"    {spectral_poly_w}")
print()


# ============================================================================
# [12] HISTORICAL SPECTRAL CENTER
# ============================================================================

print()
print("[12] HISTORICAL CENTER")
print("-" * 78)

historical_center = U + r / 2

centered_w = ex(
    spectral_poly_w.subs(z, xi + historical_center)
)

print("  center = U + r/2")
print()
print("  R(xi + U + r/2) =")
print(f"    {centered_w}")
print()


# ============================================================================
# [13] U-DISAPPEARANCE TEST
# ============================================================================

print()
print("[13] U-DISAPPEARANCE")
print("-" * 78)

contains_U = centered_w.has(U)

print(f"  centered polynomial contains U = {contains_U}")

if not contains_U:
    print("  PASS = True")
else:
    print("  PASS = False")

print()


# ============================================================================
# [14] CLOSED U,V CENTERED FORM
# ============================================================================

print()
print("[14] CLOSED U,V CENTERED FORM")
print("-" * 78)

uv_target = ex(
    (
        (r**2 + 4 * V**2 - 4 * xi**2)**2
        - 16 * r**2 * V**2
    ) / 16
)

print("  proposed:")
print("    ((r^2+4V^2-4xi^2)^2 - 16r^2V^2)/16")
print()

certificate(
    "U,V centered polynomial",
    centered_w - uv_target,
)


# ============================================================================
# [15] CONIC INVARIANT
# ============================================================================

print()
print("[15] CONIC INVARIANT")
print("-" * 78)

print("  K = U^2 - V^2")
print()

N_uv = ex(K / r**2)

print(f"  normalized N = {N_uv}")
print()


# ============================================================================
# [16] ROOT PRODUCT STRUCTURE
# ============================================================================

print()
print("[16] ROOT PRODUCT CHANNEL")
print("-" * 78)

product0 = ex(P * Q)
product1 = ex((P + 1) * (Q + 1))

print(f"  P*Q = {product0}")
print(f"  (P+1)(Q+1) = {product1}")
print()

certificate(
    "P*Q = (U^2-V^2)/r^2",
    product0 - K / r**2,
)

shifted_product_expected = ex(
    (K + 2 * U * r + r**2) / r**2
)

certificate(
    "(P+1)(Q+1) shifted-conic expression",
    product1 - shifted_product_expected,
)


# ============================================================================
# [17] PRODUCT OF THE TWO QUADRATIC CONSTANT TERMS
# ============================================================================

print()
print("[17] PRODUCT OF QUADRATIC CONSTANT TERMS")
print("-" * 78)

product_quadratics = ex(
    N_uv * (N_uv + S_modern + 1)
)

print("  N*(N+S+1) =")
print(f"    {product_quadratics}")
print()

target_product = ex(
    K * (K + 2 * U * r + r**2) / r**4
)

certificate(
    "quadratic-product invariant",
    product_quadratics - target_product,
)


# ============================================================================
# [18] CONIC TRANSLATION LAW
# ============================================================================

print()
print("[18] CONIC TRANSLATION LAW")
print("-" * 78)

K_shifted = ex(
    (U + r)**2 - V**2
)

print("  K       = U^2 - V^2")
print("  K_shift = (U+r)^2 - V^2")
print()
print(f"  K_shift-K = {ex(K_shifted - K)}")
print()

certificate(
    "translation law",
    K_shifted - K - (2 * U * r + r**2),
)


# ============================================================================
# [19] HISTORICAL INTERPRETATION
# ============================================================================

print()
print("[19] HISTORICAL INTERPRETATION")
print("-" * 78)

print(
r"""
The historical hyperbola is

    U^2 - V^2 = K.

Its factor coordinates are

    A0 = U-V
    B0 = U+V.

The KAPPA spectral bases are

    A0/r,
    B0/r,
    (A0+r)/r,
    (B0+r)/r.

Thus the +1 KAPPA shift corresponds to

    A0 -> A0+r
    B0 -> B0+r,

which is exactly

    U -> U+r
    V -> V.

After spectral centering, U disappears from the
characteristic polynomial.

The four-base spectrum therefore separates naturally into

    U^2-V^2  -> conic invariant -> N
    V^2      -> factor-gap information
    U        -> symmetric location / factor sum.

The important question is whether this centered structure
can also be transferred to the actual KAPPA sequence
F_n without first inserting P and Q.
"""
)


# ============================================================================
# [20] FINAL SYMBOLIC AUDIT
# ============================================================================

print()
print("[20] FINAL SYMBOLIC AUDIT")
print("-" * 78)

cert_raw = zero(
    chi - sp.prod(z - b for b in bases)
)

cert_factor = zero(
    chi
    -
    ((z - P) * (z - Q))
    *
    ((z - P - 1) * (z - Q - 1))
)

cert_center = zero(
    chi_centered - centered_target
)

cert_uv = zero(
    centered_w - uv_target
)

cert_conic = zero(
    product0 - K / r**2
)

cert_translation = zero(
    K_shifted - K - (2 * U * r + r**2)
)

print(f"  four-base polynomial         = {cert_raw}")
print(f"  quadratic factorization     = {cert_factor}")
print(f"  centered S/Delta polynomial = {cert_center}")
print(f"  centered U/V polynomial     = {cert_uv}")
print(f"  conic invariant channel     = {cert_conic}")
print(f"  conic translation law       = {cert_translation}")

overall = (
    cert_raw
    and cert_factor
    and cert_center
    and cert_uv
    and cert_conic
    and cert_translation
)

print()
print(f"  OVERALL EXACT AUDIT = {overall}")


# ============================================================================
# [21] NEXT RESEARCH TARGET
# ============================================================================

print()
print("[21] NEXT RESEARCH TARGET")
print("-" * 78)

print(
r"""
The exact spectral transformation is

    r*lambda =
        U-V,
        U+V,
        U-V+r,
        U+V+r.

Centering at

    z = U+r/2+xi

removes U completely.

In normalized coordinates the same statement is

    xi = z-(S+1)/2,

    chi(xi)
      =
      ((xi+1/2)^2-Delta/4)
      ((xi-1/2)^2-Delta/4).

This suggests the next experiment should not search arbitrary
features.

Instead construct the centered KAPPA sequence symbolically and
test:

    1. whether its recurrence splits into even/odd sectors;

    2. whether those sectors depend on V^2 and K=U^2-V^2
       separately;

    3. whether one sector eliminates U;

    4. whether the resulting invariant is expressible directly
       through the historical conic equation.

The key desired bridge is:

    historical conic
        -> centered spectral operator
        -> KAPPA observable
        -> N-only quantity.

That would be a genuinely new structural connection rather
than another reformulation of prime factorization.
"""
)

print()
print("=" * 78)
print("EXPERIMENT 513 FINISHED")
print("=" * 78)