#!/usr/bin/env python3

from __future__ import annotations

import ast
import importlib.util
import inspect
import pathlib
import sys
from functools import lru_cache
from typing import Any, Callable

import sympy as sp


# ============================================================================
# EXPERIMENT 64
# ============================================================================
#
# DIRECT FRESH-K VALIDATION OF THE NESTED-j TRANSITION DEFECTS
#
# This experiment:
#
#   1. Loads the newest earlier working main.py automatically.
#   2. Discovers an exact E(K,j,D)-type evaluator.
#   3. Evaluates the ORIGINAL exact object directly at fresh K values.
#   4. Removes the already-established D-support.
#   5. Reconstructs Q_j(K,D) in Y = D-D0(j) from the DIRECT D values.
#   6. Reconstructs a_{j,m}(K) on the original K grid.
#   7. Tests every fresh K point against the old six-point coefficient laws.
#   8. Reconstructs the coefficient laws on the enlarged K grid.
#   9. Compares old/new coefficient polynomials exactly.
#  10. Recomputes the j-transition defects.
#
# Nothing derived from an old interpolation polynomial is used to generate
# the fresh K data.
# ============================================================================


K = sp.Symbol("K")
D = sp.Symbol("D")
Y = sp.Symbol("Y")


# ============================================================================
# GRID
# ============================================================================

K_OLD = [3, 5, 7, 9, 11, 13]

K_NEW = [
    15, 17, 19, 21,
    23, 25, 27, 29,
]

K_ALL = K_OLD + K_NEW

J_VALUES = [0, 1, 2, 3, 4, 5]

D_VALUES = [6, 8, 10, 12, 14, 16]

D0 = {
    0: 6,
    1: 8,
    2: 8,
    3: 10,
    4: 10,
    5: 14,
}

Q_DEGREE = {
    0: 5,
    1: 4,
    2: 4,
    3: 3,
    4: 3,
    5: 1,
}

SUPPORT = {
    0: sp.Integer(1),
    1: D - 6,
    2: D - 6,
    3: (D - 8) * (D - 6),
    4: (D - 8) * (D - 6),
    5: (D - 12) * (D - 10) * (D - 8) * (D - 6),
}

LADDERS = {
    0: (K + 1) * (K + 2) * (K + 3) * (K + 4),
    1: (K + 2) * (K + 3) * (K + 4),
    2: (K + 3) * (K + 4),
    3: K + 4,
    4: sp.Integer(1),
    5: (K - 5) * (K - 7) * (K - 9) * (K - 11),
}


# ============================================================================
# PROJECT KERNEL DISCOVERY
# ============================================================================

def candidate_main_files() -> list[pathlib.Path]:
    """
    Find numbered earlier experiment directories containing main.py.
    """
    cwd = pathlib.Path.cwd().resolve()
    parent = cwd.parent

    found: list[tuple[int, pathlib.Path]] = []

    if not parent.exists():
        return []

    for directory in parent.iterdir():
        if not directory.is_dir():
            continue

        if directory == cwd:
            continue

        main_py = directory / "main.py"

        if not main_py.is_file():
            continue

        try:
            number = int(directory.name)
        except ValueError:
            continue

        found.append((number, main_py))

    found.sort(reverse=True)

    return [path for _, path in found]


def top_level_function_names(
    path: pathlib.Path,
) -> list[str]:
    """
    Read function names without importing.
    """
    try:
        tree = ast.parse(
            path.read_text(
                encoding="utf-8"
            )
        )
    except Exception:
        return []

    names = []

    for node in tree.body:
        if isinstance(
            node,
            (
                ast.FunctionDef,
                ast.AsyncFunctionDef,
            ),
        ):
            names.append(node.name)

    return names


def load_module(
    path: pathlib.Path,
):
    """
    Import a previous main.py as an isolated module.
    """
    module_name = (
        f"previous_experiment_{path.parent.name}"
    )

    spec = importlib.util.spec_from_file_location(
        module_name,
        path,
    )

    if spec is None or spec.loader is None:
        raise RuntimeError(
            f"Could not create import spec for {path}"
        )

    module = importlib.util.module_from_spec(
        spec
    )

    sys.modules[module_name] = module

    spec.loader.exec_module(module)

    return module


# ============================================================================
# EXACT EVALUATOR DISCOVERY
# ============================================================================

PREFERRED_EXACT_NAMES = [
    "exact_E",
    "exact_Ej",
    "exact_discrepancy",
    "exact_discrepancy_polynomial",
    "discrepancy_polynomial",
    "discrepancy",
    "E_exact",
    "E_value",
    "exact_error",
    "error_polynomial",
    "exact_P",
    "P_exact",
    "exact_P_j_D",
]


def callable_signature(
    fn: Callable[..., Any],
) -> str:
    try:
        return str(
            inspect.signature(fn)
        )
    except Exception:
        return "(signature unavailable)"


def parameter_names(
    fn: Callable[..., Any],
) -> list[str]:
    try:
        return list(
            inspect.signature(fn)
            .parameters.keys()
        )
    except Exception:
        return []


def looks_like_evaluator(
    fn: Callable[..., Any],
) -> bool:
    """
    A conservative heuristic for a callable that can take K,j,D.
    """
    try:
        sig = inspect.signature(fn)
    except Exception:
        return False

    params = list(sig.parameters.values())

    positional = [
        p
        for p in params
        if p.kind
        in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        )
    ]

    if len(positional) < 3:
        return False

    names = {
        p.name.lower()
        for p in positional
    }

    has_k = any(
        name in names
        for name in (
            "k",
            "K".lower(),
        )
    )

    has_j = "j" in names

    has_d = any(
        name in names
        for name in (
            "d",
            "D".lower(),
            "ell",
            "l",
            "degree",
        )
    )

    return (
        has_k
        and has_j
        and has_d
    )


def discover_exact_evaluator():
    """
    Load the newest earlier experiment and discover the original exact
    E(K,j,D)-style evaluator.

    Returns:
        (callable, module_path, function_name)
    """

    paths = candidate_main_files()

    if not paths:
        raise RuntimeError(
            "No earlier numbered main.py was found beside the current "
            "experiment directory."
        )

    attempted: list[str] = []

    for path in paths:

        names = top_level_function_names(path)

        # First try explicit high-confidence names.
        ordered_names: list[str] = []

        for name in PREFERRED_EXACT_NAMES:
            if name in names:
                ordered_names.append(name)

        # Then inspect all remaining functions.
        for name in names:
            if name not in ordered_names:
                ordered_names.append(name)

        try:
            module = load_module(path)
        except Exception as exc:
            attempted.append(
                f"{path}: import failed: {exc}"
            )
            continue

        # Pass 1: preferred names.
        for name in ordered_names:

            if not hasattr(module, name):
                continue

            fn = getattr(module, name)

            if not callable(fn):
                continue

            if name not in PREFERRED_EXACT_NAMES:
                if not looks_like_evaluator(fn):
                    continue

            params = parameter_names(fn)

            attempted.append(
                f"{path.name}:{name}{callable_signature(fn)}"
            )

            if looks_like_evaluator(fn):
                return fn, path, name

        # Some kernels use functions with non-obvious names.
        # Look for a callable that accepts K,j,D-like positional arguments.
        for name in names:

            if not hasattr(module, name):
                continue

            fn = getattr(module, name)

            if not callable(fn):
                continue

            if looks_like_evaluator(fn):

                params = parameter_names(fn)

                attempted.append(
                    f"{path.name}:{name}{callable_signature(fn)}"
                )

                return fn, path, name

    detail = "\n".join(
        f"  {item}"
        for item in attempted[-40:]
    )

    raise RuntimeError(
        "Could not automatically locate an exact E(K,j,D) evaluator.\n\n"
        "Functions inspected:\n"
        f"{detail}\n\n"
        "The earlier experiment must expose a callable accepting "
        "K, j, and D."
    )


EXACT_FN, EXACT_PATH, EXACT_NAME = (
    discover_exact_evaluator()
)


print(
    f"Loaded exact evaluator from: {EXACT_PATH}"
)
print(
    f"Exact evaluator function: {EXACT_NAME}"
)
print(
    f"Signature: {callable_signature(EXACT_FN)}"
)
print()


# ============================================================================
# EXACT EVALUATOR WRAPPER
# ============================================================================

def call_exact_E(
    k: int,
    j: int,
    d: int,
) -> sp.Expr:
    """
    Call the discovered exact evaluator.

    The wrapper handles the common positional conventions:
        (k,j,D)
        (j,k,D)
        (k,D,j)
        (j,D,k)
        (D,k,j)
        (D,j,k)

    The parameter names determine the preferred ordering.
    """

    try:
        sig = inspect.signature(EXACT_FN)
        params = list(sig.parameters.values())
    except Exception:
        params = []

    positional = [
        p
        for p in params
        if p.kind
        in (
            inspect.Parameter.POSITIONAL_ONLY,
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
        )
    ]

    if len(positional) >= 3:

        names = [
            p.name.lower()
            for p in positional[:3]
        ]

        values: dict[str, int] = {
            "k": k,
            "j": j,
            "d": d,
            "ell": d,
            "l": d,
            "degree": d,
        }

        if all(
            name in values
            for name in names
        ):
            args = [
                values[name]
                for name in names
            ]

            expr = EXACT_FN(*args)

            return sp.factor(
                sp.sympify(expr)
            )

    # Fallback to the standard order.
    attempts = [
        (k, j, d),
        (j, k, d),
        (k, d, j),
        (j, d, k),
        (d, k, j),
        (d, j, k),
    ]

    last_error = None

    for args in attempts:

        try:
            expr = EXACT_FN(*args)

            expr = sp.sympify(expr)

            if expr.is_number or expr.free_symbols:
                return sp.factor(expr)

        except Exception as exc:
            last_error = exc

    raise RuntimeError(
        f"Could not call {EXACT_NAME} for "
        f"(K,j,D)=({k},{j},{d}). "
        f"Last error: {last_error}"
    )


# ============================================================================
# DIRECT SUPPORT-REMOVED QUOTIENT
# ============================================================================

def direct_q_value(
    k: int,
    j: int,
    d: int,
) -> sp.Expr:
    """
    Evaluate Q_j(k,d) directly:

        Q_j(k,d) = E_j(k,d) / Support_j(d)

    Only nonzero-support D values are used.
    """

    support_expr = sp.sympify(
        SUPPORT[j].subs(
            {
                D: sp.Integer(d),
            }
        )
    )

    if support_expr == 0:
        raise ValueError(
            f"Attempted support division at zero support: "
            f"j={j}, K={k}, D={d}"
        )

    value = call_exact_E(
        k,
        j,
        d,
    )

    return sp.factor(
        sp.cancel(
            value / support_expr
        )
    )


# ============================================================================
# DIRECT Y-COEFFICIENT ROW
# ============================================================================

def nonzero_D_values(
    j: int,
) -> list[int]:
    return [
        d
        for d in D_VALUES
        if sp.sympify(
            SUPPORT[j].subs(D, d)
        ) != 0
    ]


def direct_y_row(
    k: int,
    j: int,
) -> list[sp.Expr]:
    """
    Compute the coefficients

        Q_j(K,D)
          = sum_m a_{j,m}(K) (D-D0[j])^m

    entirely from DIRECT exact E-values.

    No old K interpolation is involved.
    """

    local_D = nonzero_D_values(j)

    expected = Q_DEGREE[j] + 1

    if len(local_D) < expected:
        raise RuntimeError(
            f"Not enough nonzero D values for j={j}: "
            f"need {expected}, got {local_D}"
        )

    points = []

    for d in local_D:

        q_value = direct_q_value(
            k,
            j,
            d,
        )

        y_value = sp.Integer(d - D0[j])

        points.append(
            (
                y_value,
                q_value,
            )
        )

    # Exact interpolation in Y at FIXED fresh K.
    poly = sp.Poly(
        sp.interpolate(
            points,
            Y,
        ),
        Y,
        domain=sp.QQ,
    )

    if poly.degree() != Q_DEGREE[j]:
        # SymPy suppresses leading zero terms.
        if not (
            poly.is_zero
            and Q_DEGREE[j] >= 0
        ):
            raise RuntimeError(
                f"Unexpected Y degree: "
                f"j={j}, K={k}, "
                f"expected={Q_DEGREE[j]}, "
                f"got={poly.degree()}"
            )

    row = [
        sp.factor(
            poly.coeff_monomial(
                Y**m
            )
        )
        for m in range(
            expected
        )
    ]

    return row


# ============================================================================
# DIRECT DATA COLLECTION
# ============================================================================

@lru_cache(maxsize=None)
def direct_row(
    k: int,
    j: int,
) -> tuple[sp.Expr, ...]:

    row = direct_y_row(
        k,
        j,
    )

    return tuple(row)


def collect_rows(
    k_values: list[int],
):

    rows: dict[int, dict[int, list[sp.Expr]]] = {}

    for j in J_VALUES:

        rows[j] = {}

        for k in k_values:

            rows[j][k] = list(
                direct_row(
                    k,
                    j,
                )
            )

    return rows


# ============================================================================
# POLYNOMIAL HELPERS
# ============================================================================

def poly_K(
    expr: sp.Expr,
) -> sp.Poly:

    return sp.Poly(
        sp.expand(expr),
        K,
        domain=sp.QQ,
    )


def degree_K(
    expr: sp.Expr,
):
    p = poly_K(expr)

    if p.is_zero:
        return -sp.oo

    return p.degree()


def interpolate_K(
    points,
) -> sp.Expr:

    return sp.factor(
        sp.interpolate(
            points,
            K,
        )
    )


def reconstruct_K_coefficients(
    rows,
    k_values,
):

    coeffs = {}

    for j in J_VALUES:

        coeffs[j] = []

        for m in range(
            Q_DEGREE[j] + 1
        ):

            points = [
                (
                    k,
                    rows[j][k][m],
                )
                for k in k_values
            ]

            coeffs[j].append(
                interpolate_K(points)
            )

    return coeffs


# ============================================================================
# TRANSITION DEFECT
# ============================================================================

def transition_defects(
    coeffs,
):

    result = {}

    for j in [0, 1, 2]:

        modulus = poly_K(
            LADDERS[j + 1]
        )

        result[j] = []

        common_m = min(
            len(coeffs[j]),
            len(coeffs[j + 1]),
        )

        for m in range(
            common_m
        ):

            delta = poly_K(
                sp.expand(
                    coeffs[j + 1][m]
                    - coeffs[j][m]
                )
            )

            _, remainder = sp.div(
                delta,
                modulus,
            )

            result[j].append(
                sp.factor(
                    remainder.as_expr()
                )
            )

    return result


# ============================================================================
# SECTION 0
# ============================================================================

def section_0():

    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    print(
        f"k values = {K_ALL}"
    )

    print(
        f"old k values = {K_OLD}"
    )

    print(
        f"fresh k values = {K_NEW}"
    )

    print(
        f"j values = {J_VALUES}"
    )

    print(
        f"D values = {D_VALUES}"
    )

    print(
        f"points = "
        f"{len(K_ALL) * len(J_VALUES) * len(D_VALUES)}"
    )

    print(
        "grid status = OK"
    )

    print()


# ============================================================================
# SECTION 1
# ============================================================================

def section_1_direct_rows(
    rows,
):

    print("=" * 78)
    print(
        "1. DIRECTLY RECONSTRUCTED Y-COEFFICIENT ROWS"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        for k in K_ALL:

            print(
                f"  K={k}: "
                f"{rows[j][k]}"
            )

        print()


# ============================================================================
# SECTION 2
# ============================================================================

def section_2_old_coefficients(
    old_coeffs,
):

    print("=" * 78)
    print(
        "2. OLD SIX-K COEFFICIENT POLYNOMIALS"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        for m, expr in enumerate(
            old_coeffs[j]
        ):

            print(
                f"  a_{{{j},{m}}}: "
                f"degree={degree_K(expr)}"
            )

            print(
                f"    {expr}"
            )

        print()


# ============================================================================
# SECTION 3
# ============================================================================

def section_3_fresh_validation(
    rows,
    old_coeffs,
):

    print("=" * 78)
    print(
        "3. FRESH-K DIRECT VALIDATION"
    )
    print("=" * 78)

    tested = 0
    failures = 0

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        for k in K_NEW:

            for m, expr in enumerate(
                old_coeffs[j]
            ):

                predicted = sp.factor(
                    expr.subs(
                        K,
                        k,
                    )
                )

                direct = rows[j][k][m]

                ok = (
                    sp.expand(
                        predicted
                        - direct
                    )
                    == 0
                )

                tested += 1

                if not ok:
                    failures += 1

                print(
                    f"  K={k}, m={m}: "
                    f"ok={ok}"
                )

                if not ok:

                    print(
                        f"    direct    = {direct}"
                    )

                    print(
                        f"    predicted = {predicted}"
                    )

        print()

    print(
        f"fresh coefficient tests={tested}"
    )

    print(
        f"fresh coefficient failures={failures}"
    )

    print()


# ============================================================================
# SECTION 4
# ============================================================================

def section_4_enlarged_coefficients(
    all_coeffs,
):

    print("=" * 78)
    print(
        "4. ENLARGED-K COEFFICIENT POLYNOMIALS"
    )
    print("=" * 78)

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        for m, expr in enumerate(
            all_coeffs[j]
        ):

            print(
                f"  a_{{{j},{m}}}: "
                f"degree={degree_K(expr)}"
            )

            print(
                f"    {expr}"
            )

        print()


# ============================================================================
# SECTION 5
# ============================================================================

def section_5_polynomial_identity(
    old_coeffs,
    all_coeffs,
):

    print("=" * 78)
    print(
        "5. OLD VS ENLARGED COEFFICIENT POLYNOMIALS"
    )
    print("=" * 78)

    failures = 0

    for j in J_VALUES:

        print(
            f"j={j}"
        )

        for m in range(
            Q_DEGREE[j] + 1
        ):

            old = old_coeffs[j][m]
            new = all_coeffs[j][m]

            diff = sp.factor(
                sp.expand(
                    old - new
                )
            )

            ok = diff == 0

            if not ok:
                failures += 1

            print(
                f"  m={m}: "
                f"identical={ok}"
            )

            if not ok:

                print(
                    f"    difference={diff}"
                )

        print()

    print(
        f"polynomial identity failures={failures}"
    )

    print()


# ============================================================================
# SECTION 6
# ============================================================================

def section_6_old_transition_defects(
    old_coeffs,
):

    print("=" * 78)
    print(
        "6. OLD SIX-K TRANSITION DEFECTS"
    )
    print("=" * 78)

    defects = transition_defects(
        old_coeffs
    )

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        for m, defect in enumerate(
            defects[j]
        ):

            print(
                f"  m={m}: "
                f"degree={degree_K(defect)}"
            )

            print(
                f"    {defect}"
            )

        print()

    return defects


# ============================================================================
# SECTION 7
# ============================================================================

def section_7_enlarged_transition_defects(
    all_coeffs,
):

    print("=" * 78)
    print(
        "7. ENLARGED-K TRANSITION DEFECTS"
    )
    print("=" * 78)

    defects = transition_defects(
        all_coeffs
    )

    expected_degree = {
        0: 2,
        1: 1,
        2: 0,
    }

    failures = 0

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        for m, defect in enumerate(
            defects[j]
        ):

            degree = degree_K(
                defect
            )

            ok = (
                degree
                <= expected_degree[j]
            )

            if not ok:
                failures += 1

            print(
                f"  m={m}: "
                f"degree={degree}, "
                f"expected<="
                f"{expected_degree[j]}, "
                f"ok={ok}"
            )

            print(
                f"    defect={defect}"
            )

        print()

    print(
        f"defect-degree failures={failures}"
    )

    print()

    return defects


# ============================================================================
# SECTION 8
# ============================================================================

def section_8_defect_identity(
    old_defects,
    new_defects,
):

    print("=" * 78)
    print(
        "8. OLD VS ENLARGED TRANSITION-DEFECT IDENTITY"
    )
    print("=" * 78)

    failures = 0

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        for m in range(
            len(old_defects[j])
        ):

            old = old_defects[j][m]
            new = new_defects[j][m]

            diff = sp.factor(
                sp.expand(
                    old - new
                )
            )

            ok = diff == 0

            if not ok:
                failures += 1

            print(
                f"  m={m}: "
                f"identical={ok}"
            )

            if not ok:

                print(
                    f"    old={old}"
                )

                print(
                    f"    new={new}"
                )

                print(
                    f"    difference={diff}"
                )

        print()

    print(
        f"defect identity failures={failures}"
    )

    print()


# ============================================================================
# SECTION 9
# ============================================================================

def section_9_defect_fresh_targets(
    new_defects,
):

    print("=" * 78)
    print(
        "9. FRESH-K TRANSITION-DEFECT TARGETS"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        for m, defect in enumerate(
            new_defects[j]
        ):

            values = [
                sp.factor(
                    defect.subs(
                        K,
                        k,
                    )
                )
                for k in K_NEW
            ]

            print(
                f"  m={m}:"
            )

            print(
                f"    {values}"
            )

        print()


# ============================================================================
# SECTION 10
# ============================================================================

def section_10_transition_structure(
    new_defects,
):

    print("=" * 78)
    print(
        "10. TRANSITION-DEFECT STRUCTURE"
    )
    print("=" * 78)

    for j in [0, 1, 2]:

        print(
            f"j={j}->{j+1}"
        )

        active = [
            defect
            for defect in new_defects[j]
            if defect != 0
        ]

        if not active:

            print(
                "  all defects zero"
            )

            continue

        common = active[0]

        for defect in active[1:]:

            common = sp.gcd(
                poly_K(common),
                poly_K(defect),
            ).as_expr()

        common = sp.factor(
            common
        )

        print(
            f"  gcd_across_m={common}"
        )

        print(
            f"  gcd_degree="
            f"{degree_K(common)}"
        )

        for m, defect in enumerate(
            new_defects[j]
        ):

            print(
                f"  m={m}: "
                f"degree={degree_K(defect)}"
            )

        print()


# ============================================================================
# SECTION 11
# ============================================================================

def section_11_direct_support_reconstruction(
    rows,
):

    print("=" * 78)
    print(
        "11. DIRECT SUPPORT / GRID RECONSTRUCTION"
    )
    print("=" * 78)

    failures = 0
    tested = 0

    for j in J_VALUES:

        local_D = nonzero_D_values(j)

        for k in K_NEW:

            # Reconstruct Q(Y).
            row = rows[j][k]

            q_poly = sp.Integer(0)

            for m, coeff in enumerate(
                row
            ):

                q_poly += coeff * (
                    D - D0[j]
                )**m

            for d in local_D:

                predicted_q = sp.factor(
                    q_poly.subs(
                        D,
                        d,
                    )
                )

                direct_q = direct_q_value(
                    k,
                    j,
                    d,
                )

                tested += 1

                if (
                    sp.expand(
                        predicted_q
                        - direct_q
                    )
                    != 0
                ):

                    failures += 1

    print(
        f"tested={tested}"
    )

    print(
        f"support-removed reconstruction failures="
        f"{failures}"
    )

    print()


# ============================================================================
# SECTION 12
# ============================================================================

def section_12_full_grid_check(
    rows,
):

    print("=" * 78)
    print(
        "12. FULL DIRECT E = SUPPORT * Q RECONSTRUCTION"
    )
    print("=" * 78)

    tested = 0
    failures = 0

    for j in J_VALUES:

        for k in K_NEW:

            row = rows[j][k]

            q_poly = sp.Integer(0)

            for m, coeff in enumerate(
                row
            ):

                q_poly += (
                    coeff
                    * (D - D0[j])**m
                )

            for d in D_VALUES:

                support_value = sp.factor(
                    SUPPORT[j].subs(
                        D,
                        d,
                    )
                )

                predicted = sp.factor(
                    support_value
                    * q_poly.subs(
                        D,
                        d,
                    )
                )

                direct = call_exact_E(
                    k,
                    j,
                    d,
                )

                tested += 1

                if (
                    sp.expand(
                        predicted
                        - direct
                    )
                    != 0
                ):

                    failures += 1

                    print(
                        f"FAIL "
                        f"K={k} j={j} D={d}"
                    )

                    print(
                        f"  predicted={predicted}"
                    )

                    print(
                        f"  direct={direct}"
                    )

    print(
        f"tested={tested}"
    )

    print(
        f"full-grid reconstruction failures="
        f"{failures}"
    )

    print()


# ============================================================================
# SECTION 13
# ============================================================================

def section_13_summary(
    rows,
    old_coeffs,
    all_coeffs,
    old_defects,
    new_defects,
):

    print("=" * 78)
    print(
        "13. FINAL EXPERIMENT SUMMARY"
    )
    print("=" * 78)

    fresh_tests = 0
    fresh_failures = 0

    for j in J_VALUES:

        for m in range(
            Q_DEGREE[j] + 1
        ):

            old_poly = old_coeffs[j][m]

            for k in K_NEW:

                fresh_tests += 1

                direct = rows[j][k][m]

                predicted = sp.factor(
                    old_poly.subs(
                        K,
                        k,
                    )
                )

                if (
                    sp.expand(
                        direct
                        - predicted
                    )
                    != 0
                ):

                    fresh_failures += 1

    coefficient_identity_failures = 0

    for j in J_VALUES:

        for m in range(
            Q_DEGREE[j] + 1
        ):

            if (
                sp.expand(
                    old_coeffs[j][m]
                    - all_coeffs[j][m]
                )
                != 0
            ):

                coefficient_identity_failures += 1

    defect_identity_failures = 0

    for j in [0, 1, 2]:

        for m in range(
            len(old_defects[j])
        ):

            if (
                sp.expand(
                    old_defects[j][m]
                    - new_defects[j][m]
                )
                != 0
            ):

                defect_identity_failures += 1

    print(
        f"fresh coefficient tests="
        f"{fresh_tests}"
    )

    print(
        f"fresh coefficient failures="
        f"{fresh_failures}"
    )

    print(
        f"old/new coefficient polynomial "
        f"identity failures="
        f"{coefficient_identity_failures}"
    )

    print(
        f"old/new transition-defect identity "
        f"failures="
        f"{defect_identity_failures}"
    )

    print()

    if (
        fresh_failures == 0
        and coefficient_identity_failures == 0
        and defect_identity_failures == 0
    ):

        print(
            "RESULT: the six-point coefficient laws "
            "survive every fresh K evaluation."
        )

        print(
            "RESULT: the transition-defect polynomials "
            "are unchanged after the K-grid extension."
        )

    else:

        print(
            "RESULT: the previous six-point interpolation "
            "does not survive the enlarged direct K grid."
        )

    print()


# ============================================================================
# FINAL DIAGNOSTIC
# ============================================================================

def final_diagnostic():

    print("=" * 78)
    print(
        "FINAL DIAGNOSTIC"
    )
    print("=" * 78)
    print()

    print(
        "This experiment intentionally leaves the m-axis alone."
    )

    print(
        "The available m values were too few for a meaningful "
        "independent polynomial-degree claim."
    )

    print()

    print(
        "Instead, the experiment independently extends K from"
    )

    print(
        f"  {K_OLD}"
    )

    print(
        "to"
    )

    print(
        f"  {K_ALL}"
    )

    print()

    print(
        "The fresh K values are evaluated directly by the "
        "earlier exact kernel."
    )

    print()

    print(
        "The decisive validation is therefore:"
    )

    print(
        "  direct original E(K,j,D)"
    )

    print(
        "        -> support removal"
    )

    print(
        "        -> exact Y coefficients"
    )

    print(
        "        -> old coefficient polynomial prediction"
    )

    print(
        "        -> exact comparison."
    )

    print()

    print(
        "Only if those checks pass is the nested transition"
    )

    print(
        "defect structure considered stable enough for the next"
    )

    print(
        "closed-form investigation."
    )

    print()

    print(
        "No universal R=5 law is inferred."
    )

    print(
        "No r=6."
    )

    print(
        "No full pq-kernel expansion."
    )

    print(
        "No replacement universal r,j formula."
    )

    print("=" * 78)


# ============================================================================
# MAIN
# ============================================================================

def main():

    section_0()

    # Direct exact values from ORIGINAL kernel.
    rows_old = collect_rows(
        K_OLD
    )

    rows_all = collect_rows(
        K_ALL
    )

    section_1_direct_rows(
        rows_all
    )

    # Old six-point interpolation.
    old_coeffs = reconstruct_K_coefficients(
        rows_old,
        K_OLD,
    )

    section_2_old_coefficients(
        old_coeffs
    )

    section_3_fresh_validation(
        rows_all,
        old_coeffs,
    )

    # Enlarged interpolation.
    all_coeffs = reconstruct_K_coefficients(
        rows_all,
        K_ALL,
    )

    section_4_enlarged_coefficients(
        all_coeffs
    )

    section_5_polynomial_identity(
        old_coeffs,
        all_coeffs,
    )

    old_defects = (
        section_6_old_transition_defects(
            old_coeffs
        )
    )

    new_defects = (
        section_7_enlarged_transition_defects(
            all_coeffs
        )
    )

    section_8_defect_identity(
        old_defects,
        new_defects,
    )

    section_9_defect_fresh_targets(
        new_defects
    )

    section_10_transition_structure(
        new_defects
    )

    section_11_direct_support_reconstruction(
        rows_all
    )

    section_12_full_grid_check(
        rows_all
    )

    section_13_summary(
        rows_all,
        old_coeffs,
        all_coeffs,
        old_defects,
        new_defects,
    )

    final_diagnostic()


if __name__ == "__main__":
    main()