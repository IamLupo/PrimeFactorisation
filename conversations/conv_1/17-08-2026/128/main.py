#!/usr/bin/env python3

"""
========================================================================
EXPERIMENT 193
Focused support reconstruction around the (ell-3k)/4 candidate
========================================================================

Main clue from Experiment 192
-----------------------------
All 9 known nonzero anchors fit

    L = (ell - 3k)/4
    U = (ell + k)/2.

This experiment focuses on small perturbations of that geometry.

Rules
-----
1. No results.txt is read.
2. No previous experiment output is read.
3. All functions are defined locally.
4. All data are generated locally.
5. No huge unrestricted affine search is performed.
6. We explicitly test:
       - endpoint shifts by small multiples of 1/4,
       - floor/ceil conventions,
       - reflected/complemented s coordinates,
       - support widths,
       - anchor boundary positions,
       - parity compatibility.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


# ============================================================================
# BASIC TYPES
# ============================================================================

@dataclass(frozen=True)
class Anchor:
    k: int
    ell: int
    s: int
    exact: int
    label: str


@dataclass(frozen=True)
class BoundaryFormula:
    """
    Affine formula

        (a*ell + b*k + c) / 4

    with small integer coefficients.
    """
    name: str
    a: int
    b: int
    c: int

    def value(
        self,
        k: int,
        ell: int,
    ) -> Fraction:
        return Fraction(
            self.a * ell
            + self.b * k
            + self.c,
            4,
        )


@dataclass(frozen=True)
class SupportModel:
    name: str
    lower: BoundaryFormula
    upper: BoundaryFormula
    rounding: str


@dataclass(frozen=True)
class Score:
    model: SupportModel
    inside: int
    outside: int
    boundary_hits: int
    width_sum: int
    parity_compatible: int
    parity_total: int


# ============================================================================
# ANCHORS
# ============================================================================

def generate_anchors() -> list[Anchor]:
    """
    Explicitly generate the known nonzero exact anchors.

    No file is read.
    """
    return [
        Anchor(1, 10, 2, 27, "A1"),
        Anchor(3, 20, 6, 935, "A2"),
        Anchor(7, 40, 15, -1797818, "A3"),

        Anchor(3, 7, 3, -3, "B1"),
        Anchor(3, 9, 4, 9, "B2"),
        Anchor(3, 11, 5, -16, "B3"),
        Anchor(5, 11, 5, -5, "B4"),
        Anchor(5, 19, 9, -196, "B5"),
        Anchor(9, 21, 11, -84, "B6"),
    ]


# ============================================================================
# ROUNDING
# ============================================================================

def floor_fraction(
    x: Fraction,
) -> int:
    return x.numerator // x.denominator


def ceil_fraction(
    x: Fraction,
) -> int:
    return -(
        (-x.numerator)
        // x.denominator
    )


def rounded_boundary(
    x: Fraction,
    mode: str,
) -> int:
    if mode == "floor":
        return floor_fraction(x)

    if mode == "ceil":
        return ceil_fraction(x)

    if mode == "nearest":
        floor_x = floor_fraction(x)
        ceil_x = ceil_fraction(x)

        if (
            x - floor_x
            <= ceil_x - x
        ):
            return floor_x

        return ceil_x

    raise ValueError(
        f"Unknown rounding mode: {mode}"
    )


# ============================================================================
# COORDINATE TRANSFORMS
# ============================================================================

def transform_s(
    k: int,
    ell: int,
    s: int,
) -> int:
    return s


def transform_ell_minus_s(
    k: int,
    ell: int,
    s: int,
) -> int:
    return ell - s


def transform_s_minus_k(
    k: int,
    ell: int,
    s: int,
) -> int:
    return s - k


def transform_ell_k_s(
    k: int,
    ell: int,
    s: int,
) -> int:
    return ell - k - s


def transform_2s_minus_ell(
    k: int,
    ell: int,
    s: int,
) -> int:
    return 2 * s - ell


def make_transformed_anchor(
    transform_name: str,
    anchor: Anchor,
) -> int:
    if transform_name == "s":
        return transform_s(
            anchor.k,
            anchor.ell,
            anchor.s,
        )

    if transform_name == "ell-s":
        return transform_ell_minus_s(
            anchor.k,
            anchor.ell,
            anchor.s,
        )

    if transform_name == "s-k":
        return transform_s_minus_k(
            anchor.k,
            anchor.ell,
            anchor.s,
        )

    if transform_name == "ell-k-s":
        return transform_ell_k_s(
            anchor.k,
            anchor.ell,
            anchor.s,
        )

    if transform_name == "2s-ell":
        return transform_2s_minus_ell(
            anchor.k,
            anchor.ell,
            anchor.s,
        )

    raise ValueError(
        f"Unknown transform: {transform_name}"
    )


# ============================================================================
# FOCUSED BOUNDARY FAMILY
# ============================================================================

def make_lower_candidates() -> list[BoundaryFormula]:
    """
    Around

        (ell - 3k)/4.

    We allow small perturbations.
    """
    candidates = []

    for db in [-2, -1, 0, 1, 2]:
        for dc in [-2, -1, 0, 1, 2]:
            candidates.append(
                BoundaryFormula(
                    name=(
                        f"(ell{db:+d}k"
                        f"{dc:+d})/4"
                    ),
                    a=1,
                    b=-3 + db,
                    c=dc,
                )
            )

    return candidates


def make_upper_candidates() -> list[BoundaryFormula]:
    """
    Around

        (ell + k)/2
        = (2ell + 2k)/4.
    """
    candidates = []

    for db in [-2, -1, 0, 1, 2]:
        for dc in [-2, -1, 0, 1, 2]:
            candidates.append(
                BoundaryFormula(
                    name=(
                        f"(2ell{db:+d}k"
                        f"{dc:+d})/4"
                    ),
                    a=2,
                    b=2 + db,
                    c=dc,
                )
            )

    return candidates


# ============================================================================
# SUPPORT INTERVAL
# ============================================================================

def support_interval(
    model: SupportModel,
    k: int,
    ell: int,
) -> tuple[int, int]:
    lo_raw = model.lower.value(
        k,
        ell,
    )

    hi_raw = model.upper.value(
        k,
        ell,
    )

    lo = rounded_boundary(
        lo_raw,
        model.rounding,
    )

    hi = rounded_boundary(
        hi_raw,
        model.rounding,
    )

    return lo, hi


# ============================================================================
# MODEL GENERATION
# ============================================================================

def generate_models() -> list[SupportModel]:
    lower_candidates = (
        make_lower_candidates()
    )

    upper_candidates = (
        make_upper_candidates()
    )

    models = []

    for transform_name in [
        "s",
        "ell-s",
        "s-k",
        "ell-k-s",
        "2s-ell",
    ]:
        for lower in lower_candidates:
            for upper in upper_candidates:
                for rounding in [
                    "floor",
                    "ceil",
                    "nearest",
                ]:
                    models.append(
                        SupportModel(
                            name=(
                                f"{transform_name}: "
                                f"{lower.name} .. "
                                f"{upper.name} "
                                f"({rounding})"
                            ),
                            lower=lower,
                            upper=upper,
                            rounding=rounding,
                        )
                    )

    return models


# ============================================================================
# ANCHOR SCORING
# ============================================================================

def score_model(
    model: SupportModel,
    anchors: list[Anchor],
) -> Score:
    inside = 0
    outside = 0
    boundary_hits = 0
    width_sum = 0

    for anchor in anchors:
        lo, hi = support_interval(
            model,
            anchor.k,
            anchor.ell,
        )

        transformed_s = make_transformed_anchor(
            model.name.split(":")[0],
            anchor,
        )

        if (
            lo
            <= transformed_s
            <= hi
        ):
            inside += 1

            if (
                transformed_s == lo
                or transformed_s == hi
            ):
                boundary_hits += 1
        else:
            outside += 1

        width_sum += max(
            0,
            hi - lo,
        )

    parity_compatible = 0
    parity_total = 0

    # Test whether the lower and upper endpoints are integral for the
    # cases where the observed s is integral and nonzero.
    for anchor in anchors:
        lo_raw = model.lower.value(
            anchor.k,
            anchor.ell,
        )

        hi_raw = model.upper.value(
            anchor.k,
            anchor.ell,
        )

        parity_total += 2

        if lo_raw.denominator == 1:
            parity_compatible += 1

        if hi_raw.denominator == 1:
            parity_compatible += 1

    return Score(
        model=model,
        inside=inside,
        outside=outside,
        boundary_hits=boundary_hits,
        width_sum=width_sum,
        parity_compatible=parity_compatible,
        parity_total=parity_total,
    )


# ============================================================================
# RANKING
# ============================================================================

def rank_scores(
    scores: list[Score],
) -> list[Score]:
    """
    Prefer:
        1. all anchors inside,
        2. more observed endpoint hits,
        3. narrower support,
        4. more integral endpoints.
    """
    return sorted(
        scores,
        key=lambda score: (
            -score.inside,
            -score.boundary_hits,
            score.width_sum,
            -score.parity_compatible,
        ),
    )


# ============================================================================
# PRINT BEST MODELS
# ============================================================================

def print_best_models(
    ranked: list[Score],
    limit: int = 30,
) -> None:
    print()
    print("=" * 72)
    print("BEST FOCUSED SUPPORT MODELS")
    print("=" * 72)

    for index, score in enumerate(
        ranked[:limit],
        start=1,
    ):
        print()
        print(
            f"#{index:2d}"
        )

        print(
            f"  model          : "
            f"{score.model.name}"
        )

        print(
            f"  anchors inside : "
            f"{score.inside}/9"
        )

        print(
            f"  boundary hits  : "
            f"{score.boundary_hits}"
        )

        print(
            f"  width sum      : "
            f"{score.width_sum}"
        )

        print(
            f"  integral ends  : "
            f"{score.parity_compatible}/"
            f"{score.parity_total}"
        )


# ============================================================================
# EXPLICIT CANDIDATE TABLE
# ============================================================================

def print_key_models(
    anchors: list[Anchor],
) -> None:
    """
    Test the most important hand-derived support models explicitly.
    """
    formulas = [
        (
            "(ell-3k)/4 .. (ell+k)/2",
            BoundaryFormula(
                "(ell-3k)/4",
                1,
                -3,
                0,
            ),
            BoundaryFormula(
                "(ell+k)/2",
                2,
                2,
                0,
            ),
        ),
        (
            "(ell-3k+1)/4 .. (ell+k)/2",
            BoundaryFormula(
                "(ell-3k+1)/4",
                1,
                -3,
                1,
            ),
            BoundaryFormula(
                "(ell+k)/2",
                2,
                2,
                0,
            ),
        ),
        (
            "(ell-3k-1)/4 .. (ell+k)/2",
            BoundaryFormula(
                "(ell-3k-1)/4",
                1,
                -3,
                -1,
            ),
            BoundaryFormula(
                "(ell+k)/2",
                2,
                2,
                0,
            ),
        ),
        (
            "(ell-3k)/4 .. (ell+k+1)/2",
            BoundaryFormula(
                "(ell-3k)/4",
                1,
                -3,
                0,
            ),
            BoundaryFormula(
                "(ell+k+1)/2",
                2,
                2,
                1,
            ),
        ),
        (
            "(ell-3k)/4 .. (ell+k-1)/2",
            BoundaryFormula(
                "(ell-3k)/4",
                1,
                -3,
                0,
            ),
            BoundaryFormula(
                "(ell+k-1)/2",
                2,
                2,
                -1,
            ),
        ),
        (
            "ell/2-k .. ell/2+k",
            BoundaryFormula(
                "ell/2-k",
                2,
                -4,
                0,
            ),
            BoundaryFormula(
                "ell/2+k",
                2,
                4,
                0,
            ),
        ),
    ]

    print()
    print("=" * 72)
    print("KEY HAND-DERIVED MODELS")
    print("=" * 72)

    for (
        name,
        lower,
        upper,
    ) in formulas:

        for rounding in [
            "floor",
            "ceil",
            "nearest",
        ]:
            model = SupportModel(
                name="s: " + name,
                lower=lower,
                upper=upper,
                rounding=rounding,
            )

            score = score_model(
                model,
                anchors,
            )

            print(
                f"{name:38s} "
                f"{rounding:7s} "
                f"inside={score.inside}/9 "
                f"boundary={score.boundary_hits} "
                f"width={score.width_sum}"
            )


# ============================================================================
# ANCHOR POSITIONS WITH THE MAIN CANDIDATE
# ============================================================================

def print_main_candidate_positions(
    anchors: list[Anchor],
) -> None:
    """
    Show exactly where every known anchor lands under

        L=(ell-3k)/4
        U=(ell+k)/2.

    This is the main hypothesis coming out of Experiment 192.
    """
    lower = BoundaryFormula(
        "(ell-3k)/4",
        1,
        -3,
        0,
    )

    upper = BoundaryFormula(
        "(ell+k)/2",
        2,
        2,
        0,
    )

    print()
    print("=" * 72)
    print("MAIN SUPPORT HYPOTHESIS")
    print("=" * 72)

    for rounding in [
        "floor",
        "ceil",
        "nearest",
    ]:
        print()
        print(
            f"rounding = {rounding}"
        )

        model = SupportModel(
            name="main",
            lower=lower,
            upper=upper,
            rounding=rounding,
        )

        for anchor in anchors:
            lo, hi = support_interval(
                model,
                anchor.k,
                anchor.ell,
            )

            position = anchor.s - lo
            width = hi - lo

            print(
                f"  {anchor.label}: "
                f"k={anchor.k:2d} "
                f"ell={anchor.ell:2d} "
                f"s={anchor.s:3d} "
                f"support=[{lo:3d},{hi:3d}] "
                f"position={position:3d}/"
                f"{width:3d}"
            )


# ============================================================================
# FRESH SUPPORT GRID
# ============================================================================

def generate_fresh_grid() -> list[tuple[int, int, int]]:
    """
    Generate a fresh grid for examining the candidate support shape.

    This is not used as exact data; it is purely geometric.
    """
    grid = []

    for k in [
        1,
        3,
        5,
        7,
        9,
        11,
        13,
    ]:
        for ell in range(
            max(7, k),
            61,
        ):
            lo = floor_fraction(
                Fraction(
                    ell - 3 * k,
                    4,
                )
            )

            hi = floor_fraction(
                Fraction(
                    ell + k,
                    2,
                )
            )

            for s in range(
                lo,
                hi + 1,
            ):
                grid.append(
                    (
                        k,
                        ell,
                        s,
                    )
                )

    return grid


def print_width_patterns() -> None:
    """
    Print widths for fresh parameter values.
    """
    print()
    print("=" * 72)
    print("FRESH WIDTH PATTERN")
    print("=" * 72)

    for k in [
        1,
        3,
        5,
        7,
        9,
    ]:
        print()
        print(
            f"k={k}"
        )

        for ell in [
            8,
            9,
            10,
            11,
            12,
            15,
            20,
            30,
            40,
        ]:
            lo_raw = Fraction(
                ell - 3 * k,
                4,
            )

            hi_raw = Fraction(
                ell + k,
                2,
            )

            lo = floor_fraction(
                lo_raw
            )

            hi = floor_fraction(
                hi_raw
            )

            print(
                f"  ell={ell:2d}: "
                f"L={lo:3d}, "
                f"U={hi:3d}, "
                f"width={hi-lo:3d}"
            )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 72)
    print("EXPERIMENT 193")
    print("Focused support reconstruction")
    print("=" * 72)

    print()
    print(
        "All data are generated locally."
    )

    print(
        "No results.txt or previous experiment output is read."
    )

    anchors = generate_anchors()

    print()
    print(
        f"Known nonzero anchors: {len(anchors)}"
    )

    print_key_models(
        anchors
    )

    print_main_candidate_positions(
        anchors
    )

    models = generate_models()

    print()
    print(
        f"Focused models generated: "
        f"{len(models)}"
    )

    scores = [
        score_model(
            model,
            anchors,
        )
        for model in models
    ]

    ranked = rank_scores(
        scores
    )

    print_best_models(
        ranked,
        limit=30,
    )

    grid = generate_fresh_grid()

    print()
    print(
        f"Fresh geometric grid points: "
        f"{len(grid)}"
    )

    print_width_patterns()

    print()
    print("=" * 72)
    print("END EXPERIMENT 193")
    print("=" * 72)


if __name__ == "__main__":
    main()