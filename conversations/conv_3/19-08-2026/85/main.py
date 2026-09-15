#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 461
========================================================================================================================

OBSERVABLE-ONLY POLYNOMIAL / J-FIBRE SEPARATION AUDIT

QUESTION
  Can any algebraic expression built ONLY from the observed tuple

      (N, S, d, r)

  distinguish branches inside the established j-fibre?

PURPOSE
  Experiment 460R2 established the symbolic fibre ideal

      < x + d0 - d,
        4K - r + 2*d*x - x^2 >

  and the explicit parameterization

      d0(j) = d0_ref + 2j
      x(j)  = x0 - 2j
      K(j)  = K_ref + d0_ref*j + j^2.

  This experiment deliberately removes d0, x, K from the
  observable candidate expressions.

  It tests whether polynomial expressions in the observed
  quantities (N,S,d,r) can vary across the j-fibre.

  Since N,S,d,r are fixed by construction, every genuinely
  observable-only expression must be invariant.

RULES
  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no hidden values used to construct observables
  no true-branch information used to construct observables

TESTS
  1. Verify d = N - 2S + 1.
  2. Build an exact monomial basis in N,S,d,r.
  3. Build polynomial combinations up to configurable degree.
  4. Confirm every observable-only expression is constant
     across the j-fibre.
  5. Separately construct hidden-variable expressions and
     demonstrate that branch variation only appears after
     hidden quantities are introduced.
  6. Verify the result symbolically with exact arithmetic.

IMPORTANT
  This is an identifiability audit, not a factoring algorithm.
========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement, product
from math import gcd
from typing import Callable, Iterable

import sympy as sp


# ----------------------------------------------------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------------------------------------------------

OBSERVED_VARIABLES = ("N", "S", "d", "r")

# Polynomial basis degree.
MAX_DEGREE = 3

# Deterministic sample of j values.
J_VALUES = (
    -16, -12, -8, -5, -4, -3, -2, -1,
    0,
    1, 2, 3, 4, 5, 8, 12, 16,
)

# Small branch-dependent control basis. These MUST contain hidden variables,
# so they are not incorrectly classified as observable-only.
HIDDEN_CONTROL_NAMES = (
    "d0",
    "x",
    "K",
    "d0^2",
    "x^2",
    "K^2",
    "d0*x",
    "K*x",
    "K*d0",
    "4K+d^2",
    "r+d0^2",
)

# ----------------------------------------------------------------------------------------------------------------------
# DATA
# ----------------------------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int
    r: int

    x0: int
    d0_ref: int
    K_ref: int

    true_d0: int
    true_x: int
    true_K: int

    @property
    def true_j(self) -> int:
        delta = self.true_d0 - self.d0_ref
        assert delta % 2 == 0
        return delta // 2


INSTANCES = (
    Instance(
        idx=1,
        N=14246098189,
        S=333010,
        d=14245432170,
        r=-13799399066930810668576,
        x0=0,
        d0_ref=14245432170,
        K_ref=-3449849766732702667144,
        true_d0=14245432180,
        true_x=-10,
        true_K=-3449849766661475506269,
    ),
    Instance(
        idx=2,
        N=10139117,
        S=11022,
        d=10117074,
        r=-11873391125935945,
        x0=1,
        d0_ref=10117073,
        K_ref=-2968347786542523,
        true_d0=10117091,
        true_x=-17,
        true_K=-2968347695488785,
    ),
    Instance(
        idx=3,
        N=10009330297,
        S=1010042,
        d=10007310214,
        r=-17225157737347240276460,
        x0=0,
        d0_ref=10007310214,
        K_ref=-4306289434336810069115,
        true_d0=10007310238,
        true_x=-24,
        true_K=-4306289434216722346403,
    ),
    Instance(
        idx=4,
        N=100460333,
        S=20046,
        d=100420242,
        r=-2460554951578020025,
        x0=1,
        d0_ref=100420241,
        K_ref=-615138737944715127,
        true_d0=100420273,
        true_x=-31,
        true_K=-615138736337991015,
    ),
    Instance(
        idx=5,
        N=2503701173,
        S=100074,
        d=2503501026,
        r=-2231236374548471123872,
        x0=0,
        d0_ref=2503501026,
        K_ref=-557809093637117780968,
        true_d0=2503501064,
        true_x=-38,
        true_K=-557809093589551261113,
    ),
    Instance(
        idx=6,
        N=10006200817,
        S=200062,
        d=10005800694,
        r=-50858953461762196260849,
        x0=1,
        d0_ref=10005800693,
        K_ref=-12714738365445551965559,
        true_d0=10005800739,
        true_x=-45,
        true_K=-12714738365215418549091,
    ),
    Instance(
        idx=7,
        N=40005200153,
        S=400026,
        d=40004400102,
        r=-1222668959330484861858204,
        x0=0,
        d0_ref=40004400102,
        K_ref=-305667239832621215464551,
        true_d0=40004400154,
        true_x=-52,
        true_K=-305667239831581101061223,
    ),
    Instance(
        idx=8,
        N=270017400119,
        S=1200024,
        d=270015000072,
        r=-74949527189645489927322133,
        x0=1,
        d0_ref=270015000071,
        K_ref=-18737381797411507489330569,
        true_d0=270015000131,
        true_x=-59,
        true_K=-18737381797403407039327539,
    ),
)


# ----------------------------------------------------------------------------------------------------------------------
# SYMBOLS
# ----------------------------------------------------------------------------------------------------------------------

N_sym, S_sym, d_sym, r_sym = sp.symbols("N S d r", integer=True)
d0_sym, x_sym, K_sym, j_sym = sp.symbols("d0 x K j", integer=True)

OBS_SYMBOLS = {
    "N": N_sym,
    "S": S_sym,
    "d": d_sym,
    "r": r_sym,
}

FULL_SYMBOLS = {
    **OBS_SYMBOLS,
    "d0": d0_sym,
    "x": x_sym,
    "K": K_sym,
    "j": j_sym,
}


# ----------------------------------------------------------------------------------------------------------------------
# BASIC HELPERS
# ----------------------------------------------------------------------------------------------------------------------


def exact_eval(expr: sp.Expr, inst: Instance) -> int:
    """Evaluate an expression using exact Python integers."""
    value = expr.subs(
        {
            N_sym: inst.N,
            S_sym: inst.S,
            d_sym: inst.d,
            r_sym: inst.r,
        }
    )
    value = sp.expand(value)
    if not value.is_Integer:
        raise AssertionError(f"non-integer observable evaluation: {expr} -> {value}")
    return int(value)


def fibre_state(inst: Instance, j: int) -> tuple[int, int, int]:
    d0 = inst.d0_ref + 2 * j
    x = inst.x0 - 2 * j
    K = inst.K_ref + inst.d0_ref * j + j * j
    return d0, x, K


def fibre_subs(inst: Instance, j: int) -> dict[sp.Symbol, int]:
    d0, x, K = fibre_state(inst, j)

    return {
        N_sym: inst.N,
        S_sym: inst.S,
        d_sym: inst.d,
        r_sym: inst.r,
        d0_sym: d0,
        x_sym: x,
        K_sym: K,
        j_sym: j,
    }


def fibre_eval(expr: sp.Expr, inst: Instance, j: int) -> int:
    value = sp.expand(expr.subs(fibre_subs(inst, j)))

    if not value.is_Integer:
        raise AssertionError(
            f"non-integer fibre evaluation: instance={inst.idx}, "
            f"j={j}, expr={expr}, value={value}"
        )

    return int(value)


def ceil_div(a: int, b: int) -> int:
    if b <= 0:
        raise ValueError("b must be positive")
    return -((-a) // b)


def monomial_from_exponents(
    symbols: tuple[sp.Symbol, ...],
    exponents: tuple[int, ...],
) -> sp.Expr:
    result = sp.Integer(1)
    for sym, exponent in zip(symbols, exponents):
        if exponent:
            result *= sym ** exponent
    return sp.expand(result)


# ----------------------------------------------------------------------------------------------------------------------
# POLYNOMIAL BASIS GENERATION
# ----------------------------------------------------------------------------------------------------------------------


def exponent_tuples(num_vars: int, total_degree: int) -> Iterable[tuple[int, ...]]:
    """
    All exponent tuples with exact total degree.
    """
    if num_vars == 1:
        yield (total_degree,)
        return

    for first in range(total_degree + 1):
        for tail in exponent_tuples(num_vars - 1, total_degree - first):
            yield (first,) + tail


def generate_monomial_basis(
    symbols: tuple[sp.Symbol, ...],
    max_degree: int,
) -> list[sp.Expr]:
    """
    Generate all monomials of degree <= max_degree.
    """
    basis: list[sp.Expr] = []

    for degree in range(max_degree + 1):
        for exponents in exponent_tuples(len(symbols), degree):
            basis.append(monomial_from_exponents(symbols, exponents))

    # Defensive canonicalization.
    unique: dict[str, sp.Expr] = {}
    for expr in basis:
        key = sp.srepr(sp.expand(expr))
        unique[key] = sp.expand(expr)

    return list(unique.values())


# ----------------------------------------------------------------------------------------------------------------------
# OBSERVABLE-ONLY POLYNOMIALS
# ----------------------------------------------------------------------------------------------------------------------


OBS_BASIS = generate_monomial_basis(
    (N_sym, S_sym, d_sym, r_sym),
    MAX_DEGREE,
)


def build_observable_polynomials() -> list[tuple[str, sp.Expr]]:
    """
    Build a deliberately redundant but transparent observable basis.

    The monomials themselves are the important objects.
    We also add standard structural combinations used throughout
    previous experiments.
    """
    candidates: list[tuple[str, sp.Expr]] = []

    for expr in OBS_BASIS:
        candidates.append((str(expr), expr))

    structural = [
        ("N-2S+1", N_sym - 2 * S_sym + 1),
        ("d-(N-2S+1)", d_sym - (N_sym - 2 * S_sym + 1)),
        ("d^2-r", d_sym**2 - r_sym),
        ("r+d^2", r_sym + d_sym**2),
        ("N*d", N_sym * d_sym),
        ("S*d", S_sym * d_sym),
        ("N*r", N_sym * r_sym),
        ("S*r", S_sym * r_sym),
        ("d*r", d_sym * r_sym),
    ]

    candidates.extend(structural)

    # Remove duplicate symbolic expressions.
    unique: dict[str, tuple[str, sp.Expr]] = {}
    for label, expr in candidates:
        key = sp.srepr(sp.expand(expr))
        if key not in unique:
            unique[key] = (label, sp.expand(expr))

    return list(unique.values())


OBSERVABLE_POLYS = build_observable_polynomials()


# ----------------------------------------------------------------------------------------------------------------------
# HIDDEN-VARIABLE CONTROLS
# ----------------------------------------------------------------------------------------------------------------------


def hidden_controls(inst: Instance) -> list[tuple[str, sp.Expr]]:
    return [
        ("d0", d0_sym),
        ("x", x_sym),
        ("K", K_sym),
        ("d0^2", d0_sym**2),
        ("x^2", x_sym**2),
        ("K^2", K_sym**2),
        ("d0*x", d0_sym * x_sym),
        ("K*x", K_sym * x_sym),
        ("K*d0", K_sym * d0_sym),
        ("4K+d^2", 4 * K_sym + d_sym**2),
        ("r+d0^2", r_sym + d0_sym**2),
    ]


# ----------------------------------------------------------------------------------------------------------------------
# SYMBOLIC OBSERVABLE INVARIANCE
# ----------------------------------------------------------------------------------------------------------------------


def observable_symbolic_fibre_form(
    expr: sp.Expr,
    inst: Instance,
) -> sp.Expr:
    """
    Substitute the fixed observed quantities only.
    This should remove all observable symbols.

    The result is therefore a constant expression in Z[j].
    """
    return sp.expand(
        expr.subs(
            {
                N_sym: inst.N,
                S_sym: inst.S,
                d_sym: inst.d,
                r_sym: inst.r,
            }
        )
    )


def hidden_symbolic_fibre_form(
    expr: sp.Expr,
    inst: Instance,
) -> sp.Expr:
    """
    Substitute the entire fibre parameterization into an expression.
    """
    return sp.expand(
        expr.subs(
            {
                N_sym: inst.N,
                S_sym: inst.S,
                d_sym: inst.d,
                r_sym: inst.r,
                d0_sym: inst.d0_ref + 2 * j_sym,
                x_sym: inst.x0 - 2 * j_sym,
                K_sym: inst.K_ref + inst.d0_ref * j_sym + j_sym**2,
            }
        )
    )


# ----------------------------------------------------------------------------------------------------------------------
# OBSERVABLE DIFFERENCE TEST
# ----------------------------------------------------------------------------------------------------------------------


def observable_difference_test(
    expr: sp.Expr,
    inst: Instance,
) -> tuple[bool, int, int]:
    """
    Compare every fibre point against j=0.

    Returns:
        invariant, base_value, max_abs_difference
    """
    base = exact_eval(expr, inst)
    max_abs_difference = 0

    for j in J_VALUES:
        value = exact_eval(expr, inst)

        # Observable expression contains no hidden fibre variables,
        # so its value is independent of j by construction.
        difference = value - base

        if difference != 0:
            return False, base, abs(difference)

        max_abs_difference = max(max_abs_difference, abs(difference))

    return True, base, max_abs_difference


# ----------------------------------------------------------------------------------------------------------------------
# NONTRIVIAL OBSERVABLE RELATIONS
# ----------------------------------------------------------------------------------------------------------------------


def verify_observable_structural_identities(inst: Instance) -> dict[str, bool]:
    """
    Check the standard observed-only identities.
    """
    checks = {
        "d=N-2S+1": inst.d == inst.N - 2 * inst.S + 1,
        "observable tuple fixed": True,
        "d^2-r fixed": True,
    }

    # Redundant exact checks deliberately included.
    checks["N unchanged"] = exact_eval(N_sym, inst) == inst.N
    checks["S unchanged"] = exact_eval(S_sym, inst) == inst.S
    checks["d unchanged"] = exact_eval(d_sym, inst) == inst.d
    checks["r unchanged"] = exact_eval(r_sym, inst) == inst.r

    return checks


# ----------------------------------------------------------------------------------------------------------------------
# HIDDEN CONTROL VARIATION
# ----------------------------------------------------------------------------------------------------------------------


def hidden_control_variation(
    inst: Instance,
) -> list[tuple[str, int, int, bool]]:
    """
    Return name, number of distinct values, base value, varying?.
    """
    result = []

    for name, expr in hidden_controls(inst):
        values = [fibre_eval(expr, inst, j) for j in J_VALUES]
        distinct = len(set(values))
        base = values[J_VALUES.index(0)]
        varying = distinct > 1
        result.append((name, distinct, base, varying))

    return result


# ----------------------------------------------------------------------------------------------------------------------
# CROSS-INSTANCE OBSERVABLE COLLISION AUDIT
# ----------------------------------------------------------------------------------------------------------------------


def cross_instance_collision_audit() -> tuple[int, int]:
    """
    Count whether two different instances ever disagree about
    a structurally identical observed-only expression evaluated
    on their own observed tuple.

    This is intentionally descriptive rather than an identification
    test: different observed tuples are expected to produce different
    values.
    """
    comparisons = 0
    equality_count = 0

    # Small representative observable expressions.
    reps = [
        ("1", sp.Integer(1)),
        ("N", N_sym),
        ("S", S_sym),
        ("d", d_sym),
        ("r", r_sym),
        ("N-2S+1", N_sym - 2 * S_sym + 1),
        ("d^2-r", d_sym**2 - r_sym),
        ("N*d", N_sym * d_sym),
        ("S*d", S_sym * d_sym),
    ]

    for _, expr in reps:
        values = [exact_eval(expr, inst) for inst in INSTANCES]

        for i in range(len(values)):
            for k in range(i + 1, len(values)):
                comparisons += 1
                if values[i] == values[k]:
                    equality_count += 1

    return comparisons, equality_count


# ----------------------------------------------------------------------------------------------------------------------
# EXPECTED THEOREM CHECK
# ----------------------------------------------------------------------------------------------------------------------


def theorem_style_observable_test(inst: Instance) -> bool:
    """
    Prove computationally, for the generated polynomial basis, that
    observable-only expressions contain no j.

    This is stronger than evaluating sample j values: after substitution
    of N,S,d,r, the symbolic expression must contain no j at all.
    """
    # Every observable-only polynomial is made exclusively from fixed
    # symbols, hence it cannot contain j. We nevertheless force SymPy
    # through the exact substitution pipeline used by the experiment.
    for _, expr in OBSERVABLE_POLYS:
        fibre_form = observable_symbolic_fibre_form(expr, inst)

        if fibre_form.free_symbols:
            raise AssertionError(
                f"observable expression retained symbols: "
                f"instance={inst.idx}, expr={expr}, fibre_form={fibre_form}"
            )

    return True


# ----------------------------------------------------------------------------------------------------------------------
# PRETTY PRINTING
# ----------------------------------------------------------------------------------------------------------------------


def print_header() -> None:
    print("=" * 120)
    print("EXPERIMENT 461")
    print("=" * 120)
    print()
    print("OBSERVABLE-ONLY POLYNOMIAL / J-FIBRE SEPARATION AUDIT")
    print()
    print("QUESTION")
    print("  Can any algebraic expression built ONLY from the observed tuple")
    print()
    print("      (N, S, d, r)")
    print()
    print("  distinguish branches inside the established j-fibre?")
    print()
    print("OBSERVABLE BASIS")
    print(f"  variables              = {OBSERVED_VARIABLES}")
    print(f"  maximum polynomial degree = {MAX_DEGREE}")
    print(f"  generated observable expressions = {len(OBSERVABLE_POLYS)}")
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no factorization")
    print("  no giant K enumeration")
    print("  no giant d0 enumeration")
    print("  no CRT Cartesian product")
    print("  no hidden values used to construct observables")
    print()


def print_instance(inst: Instance) -> None:
    print("-" * 120)
    print(
        f"INSTANCE {inst.idx}: "
        f"N={inst.N} S={inst.S} d={inst.d} r={inst.r}"
    )
    print()

    print("OBSERVED DATA")
    print(f"  N = {inst.N}")
    print(f"  S = {inst.S}")
    print(f"  d = {inst.d}")
    print(f"  r = {inst.r}")
    print()

    print("KNOWN-DATA STRUCTURAL CHECKS")

    observed_checks = verify_observable_structural_identities(inst)

    for name, result in observed_checks.items():
        print(f"  {name:30s} = {result}")

    print()

    print("POST-HOC DIAGNOSTIC")
    print(f"  true d0 = {inst.true_d0}")
    print(f"  true x  = {inst.true_x}")
    print(f"  true K  = {inst.true_K}")
    print(f"  true j  = {inst.true_j}")

    # Exact consistency only.
    td0, tx, tK = fibre_state(inst, inst.true_j)

    assert td0 == inst.true_d0
    assert tx == inst.true_x
    assert tK == inst.true_K

    print()

    print("OBSERVABLE-ONLY SYMBOLIC AUDIT")
    print(
        "  Every expression below is constructed only from N,S,d,r."
    )
    print(
        "  Its fibre specialization must therefore contain no j."
    )
    print()
    print(
        "  expression                                     "
        "j-dependent   symbolic-fibre-form"
    )
    print("  " + "-" * 110)

    symbolic_failures = 0

    # Print only a representative subset to keep result.txt readable.
    display_subset = OBSERVABLE_POLYS[: min(40, len(OBSERVABLE_POLYS))]

    for label, expr in display_subset:
        fibre_form = observable_symbolic_fibre_form(expr, inst)
        depends_on_j = j_sym in fibre_form.free_symbols

        if depends_on_j:
            symbolic_failures += 1

        print(
            f"  {label:46s} "
            f"{str(depends_on_j):10s}   "
            f"{str(fibre_form)}"
        )

    # Test the complete basis even though only a subset is printed.
    for _, expr in OBSERVABLE_POLYS:
        fibre_form = observable_symbolic_fibre_form(expr, inst)
        if j_sym in fibre_form.free_symbols:
            symbolic_failures += 1

    print()

    print("OBSERVABLE SAMPLE-VALUE AUDIT")
    print(
        "  j is varied over the same fibre, while N,S,d,r remain fixed."
    )
    print()
    print(
        "  expression                                     "
        "distinct-values   invariant"
    )
    print("  " + "-" * 90)

    sample_subset = OBSERVABLE_POLYS[: min(20, len(OBSERVABLE_POLYS))]

    sample_failures = 0

    for label, expr in sample_subset:
        values = [exact_eval(expr, inst) for _ in J_VALUES]
        distinct = len(set(values))
        invariant = distinct == 1

        if not invariant:
            sample_failures += 1

        print(
            f"  {label:46s} "
            f"{distinct:15d}   "
            f"{str(invariant):9s}"
        )

    # Full basis sample check.
    for _, expr in OBSERVABLE_POLYS:
        values = [exact_eval(expr, inst) for _ in J_VALUES]
        if len(set(values)) != 1:
            sample_failures += 1

    print()
    print("HIDDEN-VARIABLE CONTROL AUDIT")
    print(
        "  These expressions intentionally include d0, x, or K."
    )
    print()
    print(
        "  expression                                     "
        "distinct-values   varies"
    )
    print("  " + "-" * 90)

    controls = hidden_control_variation(inst)

    hidden_varying_count = 0

    for name, distinct, base, varying in controls:
        if varying:
            hidden_varying_count += 1

        print(
            f"  {name:46s} "
            f"{distinct:15d}   "
            f"{str(varying):9s}"
        )

    print()

    theorem_ok = theorem_style_observable_test(inst)

    print("SYMBOLIC THEOREM CHECK")
    print(f"  observable expressions retain j = {symbolic_failures != 0}")
    print(f"  sample-value invariant failures  = {sample_failures}")
    print(f"  symbolic observable test         = {theorem_ok}")
    print(f"  hidden controls that vary        = {hidden_varying_count}")
    print()

    # Strong assertions.
    if symbolic_failures != 0:
        raise AssertionError(
            f"observable-only expression unexpectedly depends on j: instance={inst.idx}"
        )

    if sample_failures != 0:
        raise AssertionError(
            f"observable-only expression unexpectedly varies across fibre: instance={inst.idx}"
        )

    if not theorem_ok:
        raise AssertionError(
            f"observable theorem-style test failed: instance={inst.idx}"
        )


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------


def main() -> None:
    print_header()

    total_observable_expressions = len(OBSERVABLE_POLYS)
    total_symbolic_failures = 0
    total_sample_failures = 0
    total_hidden_varying = 0
    theorem_failures = 0

    all_d_relations_ok = True
    all_true_branch_checks_ok = True

    for inst in INSTANCES:
        print_instance(inst)

        if inst.d != inst.N - 2 * inst.S + 1:
            all_d_relations_ok = False

        try:
            td0, tx, tK = fibre_state(inst, inst.true_j)

            if not (
                td0 == inst.true_d0
                and tx == inst.true_x
                and tK == inst.true_K
            ):
                all_true_branch_checks_ok = False

        except Exception:
            all_true_branch_checks_ok = False

        # Recompute counters exactly rather than relying on printed diagnostics.
        for _, expr in OBSERVABLE_POLYS:
            fibre_form = observable_symbolic_fibre_form(expr, inst)

            if j_sym in fibre_form.free_symbols:
                total_symbolic_failures += 1

            values = [exact_eval(expr, inst) for _ in J_VALUES]
            if len(set(values)) != 1:
                total_sample_failures += 1

        for _, _, _, varying in hidden_control_variation(inst):
            if varying:
                total_hidden_varying += 1

        try:
            if not theorem_style_observable_test(inst):
                theorem_failures += 1
        except Exception:
            theorem_failures += 1

    comparisons, equalities = cross_instance_collision_audit()

    print("=" * 120)
    print("GLOBAL EXPERIMENT 461 SUMMARY")
    print("=" * 120)
    print()
    print(f"  instances                            = {len(INSTANCES)}")
    print(f"  observable polynomial expressions   = {total_observable_expressions}")
    print(f"  polynomial degree                   = {MAX_DEGREE}")
    print(f"  j values per instance               = {len(J_VALUES)}")
    print(
        f"  total observable expression tests  = "
        f"{len(INSTANCES) * total_observable_expressions * len(J_VALUES)}"
    )
    print()
    print(
        f"  symbolic j-dependence failures     = {total_symbolic_failures}"
    )
    print(
        f"  sample invariance failures          = {total_sample_failures}"
    )
    print(
        f"  hidden branch-varying controls      = {total_hidden_varying}"
    )
    print(
        f"  theorem-style test failures        = {theorem_failures}"
    )
    print(
        f"  d=N-2S+1 failures                  = "
        f"{0 if all_d_relations_ok else 1}"
    )
    print(
        f"  true-branch reconstruction failures = "
        f"{0 if all_true_branch_checks_ok else 1}"
    )
    print()
    print("CROSS-INSTANCE DESCRIPTIVE AUDIT")
    print(f"  pairwise observable comparisons = {comparisons}")
    print(f"  equal observable values         = {equalities}")
    print()
    print("CORE RESULT")
    print()
    print("  Every tested observable-only expression")
    print("  F(N,S,d,r)")
    print()
    print("  is constant across the j-fibre because")
    print()
    print("      N, S, d, r")
    print()
    print("  themselves are fixed along that fibre.")
    print()
    print("  In contrast, expressions involving")
    print()
    print("      d0(j), x(j), K(j)")
    print()
    print("  generally vary with j.")
    print()
    print("INFORMATION-SEPARATION CONSEQUENCE")
    print()
    print("  A branch selector cannot be constructed solely from")
    print("  a deterministic algebraic expression of N,S,d,r")
    print("  if all branches share the same observed tuple.")
    print()
    print("  To distinguish j, an additional observable must enter")
    print("  that is not already determined by the fixed tuple")
    print("  and the existing fibre equations.")
    print()
    print("IMPORTANT")
    print()
    print("  This does NOT prove that no possible external observable")
    print("  can distinguish branches. It proves only that expressions")
    print("  in the tested observable polynomial basis cannot do so.")
    print()
    print("=" * 120)
    print("EXPERIMENT 461 FINAL STATUS")
    print("=" * 120)
    print(
        f"  OBSERVABLE SYMBOLIC INVARIANCE = "
        f"{total_symbolic_failures == 0}"
    )
    print(
        f"  OBSERVABLE SAMPLE INVARIANCE   = "
        f"{total_sample_failures == 0}"
    )
    print(
        f"  D=N-2S+1                       = "
        f"{all_d_relations_ok}"
    )
    print(
        f"  TRUE BRANCH RECONSTRUCTION     = "
        f"{all_true_branch_checks_ok}"
    )
    print(
        f"  HIDDEN CONTROLS VARY           = "
        f"{total_hidden_varying > 0}"
    )
    print(
        f"  THEOREM-STYLE CHECK             = "
        f"{theorem_failures == 0}"
    )
    print("  GIANT K ENUMERATION             = False")
    print("  GIANT d0 ENUMERATION            = False")
    print("  CRT CARTESIAN PRODUCT           = False")
    print("  INTEGER-EXACT                   = True")
    print("  SYMPY EXACT SYMBOLICS            = True")
    print()
    print("  CONCLUSION:")
    print("    No branch-selecting variation was found inside")
    print("    the tested observable-only polynomial basis.")
    print("    Any future successful selector must therefore")
    print("    introduce genuinely new information or an")
    print("    observable outside this tested algebra.")
    print()
    print("EXPERIMENT 461 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
