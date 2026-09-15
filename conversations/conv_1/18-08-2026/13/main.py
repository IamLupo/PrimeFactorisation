from __future__ import annotations

import ast
import importlib.util
import pathlib
import sys
from functools import lru_cache
from typing import Callable, Optional

import sympy as sp

# ============================================================
# EXPERIMENT 225
# COMPLETE EXACT L3 BOUNDARY DERIVATION
# ============================================================
#
# This version is self-contained with respect to the PROJECT
# KERNEL: it automatically locates the newest earlier main.py
# beside the current experiment and reuses exact_F() or
# exact_Q_pq(). It does NOT require a manual placeholder.
#
# All arithmetic is exact over QQ.
# ============================================================

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")
Lsym = sp.symbols("L", integer=True, nonnegative=True)


# ============================================================
# PROJECT KERNEL LOADER
# ============================================================

def _candidate_main_files() -> list[pathlib.Path]:
    """Find earlier experiment main.py files in sibling directories."""
    cwd = pathlib.Path.cwd().resolve()
    parent = cwd.parent
    candidates: list[tuple[int, pathlib.Path]] = []

    for path in parent.iterdir():
        if not path.is_dir() or path == cwd:
            continue
        main = path / "main.py"
        if not main.is_file():
            continue
        try:
            number = int(path.name)
        except ValueError:
            continue
        candidates.append((number, main))

    candidates.sort(reverse=True)
    return [path for _, path in candidates]


def _defines_function(path: pathlib.Path, names: set[str]) -> bool:
    """Check source AST without importing a candidate module."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    return any(
        isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name in names
        for node in tree.body
    )


def _load_previous_kernel() -> tuple[Callable[[int, int], sp.Expr], pathlib.Path, str]:
    """Load exact_F or exact_Q_pq from newest earlier experiment."""
    candidates = _candidate_main_files()
    names = {"exact_F", "exact_Q_pq"}

    for path in candidates:
        if not _defines_function(path, names):
            continue

        module_name = f"prior_experiment_{path.parent.name}"
        spec = importlib.util.spec_from_file_location(module_name, path)
        if spec is None or spec.loader is None:
            continue

        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module

        try:
            spec.loader.exec_module(module)
        except Exception as exc:
            print(
                f"WARNING: could not import prior kernel from {path}: {exc}",
                file=sys.stderr,
            )
            continue

        if hasattr(module, "exact_F"):
            fn = getattr(module, "exact_F")
            if callable(fn):
                return fn, path, "exact_F"

        if hasattr(module, "exact_Q_pq"):
            fn = getattr(module, "exact_Q_pq")
            if callable(fn):
                return fn, path, "exact_Q_pq"

    raise RuntimeError(
        "Could not locate an earlier main.py defining exact_F() or exact_Q_pq().\n"
        "Run this experiment from the numbered experiment directory whose parent\n"
        "contains the earlier working experiments, or copy the working kernel\n"
        "function into this file."
    )


_KERNEL_FN, _KERNEL_PATH, _KERNEL_NAME = _load_previous_kernel()

print("Loaded exact kernel from:")
print(f"  {_KERNEL_PATH}")
print(f"  function = {_KERNEL_NAME}")
print()


def exact_F(k: int, ell: int) -> sp.Expr:
    """Stable wrapper around the project's exact kernel."""
    expr = _KERNEL_FN(k, ell)
    expr = sp.expand(sp.sympify(expr))

    # The previous experiment chain is symmetric in p,q.
    if sp.expand(expr.xreplace({p: q, q: p}) - expr) != 0:
        # Some project implementations may use their own symbols.
        # Rebuild the expression in our p,q symbols by substitution.
        free = sorted(expr.free_symbols, key=lambda s: s.name)
        if len(free) >= 2:
            sub = {free[0]: p, free[1]: q}
            expr = sp.expand(expr.xreplace(sub))

    return expr


# ============================================================
# EXACT SYMMETRIC REDUCTION
# ============================================================

@lru_cache(maxsize=None)
def exact_G(k: int, ell: int) -> sp.Expr:
    """Return G(N,X), where N=pq and X=p+q+1."""
    F = sp.expand(exact_F(k, ell))
    S = X - 1
    poly = sp.Poly(F, p, q)
    max_deg = poly.total_degree()

    power_sum: dict[int, sp.Expr] = {
        0: sp.Integer(2),
        1: S,
    }
    for n in range(2, max_deg + 1):
        power_sum[n] = sp.expand(
            S * power_sum[n - 1] - N * power_sum[n - 2]
        )

    result = sp.Integer(0)
    visited: set[tuple[int, int]] = set()

    for (a, b), coeff in poly.terms():
        if (a, b) in visited:
            continue

        if a == b:
            result += coeff * N**a
            visited.add((a, b))
            continue

        c1 = poly.coeff_monomial(p**a * q**b)
        c2 = poly.coeff_monomial(p**b * q**a)
        if sp.expand(c1 - c2) != 0:
            raise ValueError(
                f"Kernel is not symmetric at monomials ({a},{b}) and ({b},{a})."
            )

        hi, lo = max(a, b), min(a, b)
        d = hi - lo
        result += c1 * N**lo * power_sum[d]

        visited.add((a, b))
        visited.add((b, a))

    return sp.Poly(sp.expand(result), N, X).as_expr()


# ============================================================
# HOMOGENEOUS EXTRACTION
# ============================================================

def homogeneous_component(expr: sp.Expr, degree: int) -> sp.Expr:
    poly = sp.Poly(sp.expand(expr), N, X)
    out = sp.Integer(0)
    for (a, b), coeff in poly.terms():
        if a + b == degree:
            out += coeff * N**a * X**b
    return sp.expand(out)


# ============================================================
# TOP LAYER
# ============================================================

def G_top(k: int, ell: int) -> sp.Expr:
    return sp.expand(
        -X**(ell - k) * ((X + N)**k - N**k)
    )


# ============================================================
# L1 CLOSED FORM
# ============================================================

def L1_coeff(k: int, ell: int, a: int) -> sp.Expr:
    if a < k:
        return sp.expand(
            ell * sp.binomial(k + 1, a)
            - k * sp.binomial(k, a - 1)
        )
    if a == k:
        return sp.expand(ell * (k + 1) - k**2 + k)
    return sp.Integer(0)


def L1_closed(k: int, ell: int) -> sp.Expr:
    out = sp.Integer(0)
    for a in range(k + 1):
        out += L1_coeff(k, ell, a) * N**a * X**(ell - 1 - a)
    return sp.expand(out)


# ============================================================
# L2 CLOSED FORM — ESTABLISHED IN EXPERIMENT 217
# ============================================================

def L2_interior(k: int, ell: int, a: int) -> sp.Expr:
    Cka = sp.binomial(k + 2, a)
    return sp.expand(
        -sp.Rational(1, 2) * Cka * ell**2
        + Cka
        * (2 * (k + 1) * a + k + 2)
        * ell
        / (2 * (k + 2))
        - k * a * (a + 1) * Cka
        / (2 * (k + 2))
    )


def L2_coeff(k: int, ell: int, a: int) -> sp.Expr:
    if 0 <= a < k:
        return L2_interior(k, ell, a)

    if a == k:
        return sp.expand(
            L2_interior(k, ell, a)
            + sp.binomial(k + 2, 3)
        )

    if a == k + 1:
        core = (
            ell**2 * k
            + 2 * ell**2
            - 2 * ell * k**2
            - 3 * ell * k
            - 2 * ell
            + k**3
            + k**2
            - 2 * k
        )
        return sp.expand(-core / 2)

    return sp.Integer(0)


def L2_closed(k: int, ell: int) -> sp.Expr:
    out = sp.Integer(0)
    for a in range(k + 2):
        out += L2_coeff(k, ell, a) * N**a * X**(ell - 2 - a)
    return sp.expand(out)


# ============================================================
# EXACT L3
# ============================================================

@lru_cache(maxsize=None)
def exact_L3(k: int, ell: int) -> sp.Expr:
    G = exact_G(k, ell)
    residual = sp.expand(
        G
        - G_top(k, ell)
        - L1_closed(k, ell)
        - L2_closed(k, ell)
    )
    return homogeneous_component(residual, ell - 3)


def coeff_L3(L3: sp.Expr, a: int, ell: int) -> sp.Expr:
    b = ell - 3 - a
    if b < 0:
        return sp.Integer(0)
    return sp.Poly(L3, N, X).coeff_monomial(N**a * X**b)


# ============================================================
# ESTABLISHED INTERIOR L3 CUBIC
# ============================================================

def L3_interior(k: int, ell: int, a: int) -> sp.Expr:
    Cka = sp.binomial(k + 3, a)

    A = Cka / 6
    B = -Cka * (a * (k + 2) + k + 3) / (2 * (k + 3))
    C = Cka * (
        sp.Rational(1, 3)
        + a * ((k + 1) * a + 2 * k + 3)
        / (2 * (k + 3))
    )
    D = -k * a * (a + 1) * (a + 2) * Cka / (6 * (k + 3))

    # Cancel rationals before returning.
    return sp.factor(sp.cancel(
        A * ell**3 + B * ell**2 + C * ell + D
    ))


# ============================================================
# 1. SUPPORT AUDIT
# ============================================================

def support_audit() -> None:
    print("=" * 78)
    print("EXPERIMENT 225")
    print("COMPLETE EXACT L3 BOUNDARY DERIVATION")
    print("=" * 78)
    print()
    print("=" * 78)
    print("1. L3 SUPPORT AUDIT")
    print("=" * 78)

    tests = [
        (3, 7), (3, 9), (3, 11), (3, 13),
        (5, 11), (5, 13), (5, 15), (5, 17),
        (7, 15), (7, 17),
        (9, 21), (9, 23),
        (11, 23), (11, 25),
    ]

    failures = 0
    for k, ell in tests:
        L3 = exact_L3(k, ell)
        poly = sp.Poly(L3, N, X)
        support = sorted({a for (a, _), c in poly.terms() if c != 0})
        expected = list(range(k + 2))
        ok = support == expected
        print(
            f"k={k:2d} ell={ell:2d} support={support} "
            f"expected={expected} {'PASS' if ok else 'FAIL'}"
        )
        failures += int(not ok)

    print()
    print(f"support failures = {failures}")


# ============================================================
# 2. INTERIOR VALIDATION
# ============================================================

def interior_audit() -> None:
    print()
    print("=" * 78)
    print("2. EXACT INTERIOR L3 VALIDATION")
    print("=" * 78)

    tests = [
        (3, 7), (3, 9), (3, 11), (3, 13), (3, 15),
        (5, 11), (5, 13), (5, 15), (5, 17), (5, 19),
        (7, 15), (7, 17), (7, 19),
        (9, 21), (9, 23),
        (11, 23), (11, 25),
    ]

    failures = 0
    total = 0

    for k, ell in tests:
        L3 = exact_L3(k, ell)
        local_failures = []

        for a in range(k):
            total += 1
            actual = sp.cancel(coeff_L3(L3, a, ell))
            predicted = sp.cancel(L3_interior(k, ell, a))
            diff = sp.cancel(actual - predicted)
            if diff != 0:
                local_failures.append((a, actual, predicted, diff))

        if local_failures:
            failures += len(local_failures)
            for a, actual, predicted, diff in local_failures:
                print(f"FAIL k={k} ell={ell} a={a}")
                print(f"  actual    = {actual}")
                print(f"  predicted = {predicted}")
                print(f"  difference= {diff}")

        print(
            f"k={k:2d} ell={ell:2d} tested={k} "
            f"failures={len(local_failures)} "
            f"{'PASS' if not local_failures else 'FAIL'}"
        )

    print()
    print(f"tested interior coefficients = {total}")
    print(f"interior formula failures = {failures}")


# ============================================================
# 3. a=k CORRECTION
# ============================================================

def boundary_k_audit() -> None:
    print()
    print("=" * 78)
    print("3. a=k EXACT BOUNDARY CORRECTION")
    print("=" * 78)

    tests = [
        (3, [7, 9, 11, 13, 15]),
        (5, [11, 13, 15, 17, 19]),
        (7, [15, 17, 19]),
        (9, [21, 23, 25]),
        (11, [23, 25]),
        (13, [25, 27]),
    ]

    failures = 0

    for k, ells in tests:
        expected = sp.binomial(k + 2, 3)
        print()
        print(f"k={k} expected correction={expected}")

        for ell in ells:
            L3 = exact_L3(k, ell)
            actual = sp.cancel(coeff_L3(L3, k, ell))
            naive = sp.cancel(L3_interior(k, ell, k))
            correction = sp.cancel(actual - naive)
            ok = correction == expected
            print(
                f"  ell={ell:2d} actual={str(actual):>8} "
                f"correction={str(correction):>8} "
                f"expected={str(expected):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )
            failures += int(not ok)

    print()
    print(f"a=k correction failures = {failures}")


# ============================================================
# 4. a=k+1 DIRECT EXTRACTION
# ============================================================

def boundary_k1_extraction() -> None:
    print()
    print("=" * 78)
    print("4. a=k+1 DIRECT BOUNDARY EXTRACTION")
    print("=" * 78)

    tests = [
        (3, [7, 9, 11, 13, 15]),
        (5, [11, 13, 15, 17, 19]),
        (7, [15, 17, 19]),
        (9, [21, 23, 25]),
        (11, [23, 25]),
        (13, [25, 27]),
    ]

    for k, ells in tests:
        print()
        print(f"k={k}")
        points = []

        for ell in ells:
            L3 = exact_L3(k, ell)
            value = sp.cancel(coeff_L3(L3, k + 1, ell))
            points.append((ell, value))
            print(f"  ell={ell:2d} D_(k+1)={value}")

        if len(points) >= 4:
            poly = sp.factor(sp.interpolate(points, Lsym))
            print(f"  exact cubic reconstruction = {poly}")
        else:
            poly = sp.factor(sp.interpolate(points, Lsym))
            print(f"  diagnostic interpolation = {poly}")


# ============================================================
# 5. a=k+1 RELATION TO INTERIOR FORMULA
# ============================================================

def boundary_k1_correction() -> None:
    print()
    print("=" * 78)
    print("5. a=k+1 BOUNDARY CORRECTION")
    print("=" * 78)

    tests = [
        (3, [7, 9, 11, 13, 15]),
        (5, [11, 13, 15, 17, 19]),
        (7, [15, 17, 19]),
        (9, [21, 23, 25]),
        (11, [23, 25]),
    ]

    for k, ells in tests:
        print()
        print(f"k={k}")
        for ell in ells:
            L3 = exact_L3(k, ell)
            actual = sp.cancel(coeff_L3(L3, k + 1, ell))
            extrapolated = sp.cancel(L3_interior(k, ell, k + 1))
            correction = sp.cancel(actual - extrapolated)
            print(
                f"  ell={ell:2d} actual={str(actual):>8} "
                f"interior-extrap={str(extrapolated):>8} "
                f"correction={str(correction):>8}"
            )


# ============================================================
# 6. ZERO TAIL
# ============================================================

def zero_tail_audit() -> None:
    print()
    print("=" * 78)
    print("6. ZERO-TAIL AUDIT: a >= k+2")
    print("=" * 78)

    tests = [
        (3, 9),
        (5, 13),
        (7, 17),
        (9, 23),
        (11, 25),
        (13, 27),
    ]

    failures = 0

    for k, ell in tests:
        L3 = exact_L3(k, ell)
        bad = []

        for a in range(k + 2, max(k + 2, ell - 2)):
            value = sp.cancel(coeff_L3(L3, a, ell))
            if value != 0:
                bad.append((a, value))

        print(
            f"k={k:2d} ell={ell:2d} "
            f"tail={[(a, str(v)) for a, v in bad]}"
        )
        failures += len(bad)

    print()
    print(f"nonzero-tail failures = {failures}")


# ============================================================
# 7. THIRD-DIFFERENCE LAW
# ============================================================

def cubic_difference_audit() -> None:
    print()
    print("=" * 78)
    print("7. THIRD-DIFFERENCE LAW")
    print("=" * 78)

    cases = [
        (3, [7, 9, 11, 13]),
        (5, [11, 13, 15, 17]),
        (7, [15, 17, 19, 21]),
        (9, [21, 23, 25]),
        (11, [23, 25, 27]),
    ]

    failures = 0

    for k, ells in cases:
        for a in range(k):
            vals = [
                sp.cancel(
                    coeff_L3(exact_L3(k, ell), a, ell)
                )
                for ell in ells
            ]

            d1 = [vals[i + 1] - vals[i] for i in range(len(vals) - 1)]
            d2 = [d1[i + 1] - d1[i] for i in range(len(d1) - 1)]
            d3 = [d2[i + 1] - d2[i] for i in range(len(d2) - 1)]

            expected = 8 * sp.binomial(k + 3, a)
            ok = all(v == expected for v in d3)

            if not ok:
                failures += 1
                print(
                    f"FAIL k={k} a={a} "
                    f"Delta3={[str(v) for v in d3]} "
                    f"expected={expected}"
                )

    print(f"third-difference failures = {failures}")


# ============================================================
# 8. RECONSTRUCTION USING ACTUAL a=k+1 ROW
# ============================================================

def reconstruction_audit() -> None:
    print()
    print("=" * 78)
    print("8. COMPLETE L3 RECONSTRUCTION FROM DERIVED SUPPORT")
    print("=" * 78)

    tests = [
        (3, 7), (3, 9), (3, 11), (3, 13),
        (5, 11), (5, 13), (5, 15), (5, 17),
        (7, 15), (7, 17),
        (9, 21), (9, 23),
        (11, 23), (11, 25),
    ]

    failures = 0

    for k, ell in tests:
        exact = exact_L3(k, ell)
        reconstructed = sp.Integer(0)

        for a in range(k):
            reconstructed += (
                L3_interior(k, ell, a)
                * N**a
                * X**(ell - 3 - a)
            )

        reconstructed += (
            (
                L3_interior(k, ell, k)
                + sp.binomial(k + 2, 3)
            )
            * N**k
            * X**(ell - 3 - k)
        )

        boundary = coeff_L3(exact, k + 1, ell)
        reconstructed += (
            boundary
            * N**(k + 1)
            * X**(ell - 4 - k)
        )

        residual = sp.Poly(
            sp.expand(exact - reconstructed), N, X
        )
        bad = [
            (mon, coeff)
            for mon, coeff in residual.terms()
            if coeff != 0
        ]

        ok = not bad
        print(
            f"k={k:2d} ell={ell:2d} "
            f"reconstruction={'PASS' if ok else 'FAIL'}"
        )
        if bad:
            print(
                "  residual = "
                f"{[(m, str(c)) for m, c in bad]}"
            )
            failures += 1

    print()
    print(f"reconstruction failures = {failures}")


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    support_audit()
    interior_audit()
    boundary_k_audit()
    boundary_k1_extraction()
    boundary_k1_correction()
    zero_tail_audit()
    cubic_difference_audit()
    reconstruction_audit()

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print(
        """
Experiment 225 is now wired to the project's existing exact kernel.

It tests:
  * exact L3 support;
  * the established interior cubic law;
  * the exact a=k correction + C(k+2,3);
  * the a=k+1 boundary row;
  * the zero tail a>=k+2;
  * and complete L3 reconstruction.

No previous experiment output is read as data.
The previous exact kernel implementation is reused only as the
mathematical definition of F_{k,ell}(p,q).

If the a=k+1 row simplifies and the tail is zero, L3 is closed.
Do not start L4 before the complete L3 support is verified.
"""
    )


if __name__ == "__main__":
    main()
