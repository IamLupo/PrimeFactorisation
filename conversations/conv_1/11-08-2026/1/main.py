#!/usr/bin/env python3

"""
KAPPA STRUCTURE EXPLORER III
============================

Two investigations:

A) LOCAL LIFTING
----------------
For

    f(p) = p^2 - p + 1

study

    v_q(f(p))

for q, q^2, q^3, ...

We test whether valuation is completely determined by
p mod q^k and whether the roots lift uniquely.

B) GLOBAL SYMMETRIC STRUCTURE
-----------------------------
For a prime set

    P = {p1,...,px}

we have

    kappa = product f(pi)

where

    f(x) = x^2 - x + 1.

Expand:

    product (pi^2 - pi + 1)

as a symmetric polynomial in the pi.

We then compare kappa against elementary symmetric
polynomials e1,e2,...,ex.

For x=2:

    e1 = p+q
    e2 = pq=n

For x=3:

    e1 = p+q+r
    e2 = pq+pr+qr
    e3 = pqr=n

For x=4:

    e1,e2,e3,e4=n.

The goal is to discover the minimum symmetric information
needed to determine kappa.

This is NOT a factoring algorithm. It is structural
experimentation.
"""

from itertools import combinations
from math import prod
from collections import defaultdict
import csv
import json


# ============================================================
# CONFIGURATION
# ============================================================

MAX_PRIME = 200
MAX_BODIES = 4

Q_LIMIT = 50

MAX_LIFT = 5

OUTPUT_LOCAL = "kappa_lifting.csv"
OUTPUT_GLOBAL = "kappa_symmetric.csv"
OUTPUT_SUMMARY = "kappa_structure3_summary.json"


# ============================================================
# PRIME GENERATION
# ============================================================

def is_prime(n):

    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3

    while d * d <= n:

        if n % d == 0:
            return False

        d += 2

    return True


def primes_up_to(limit):

    return [
        n for n in range(2, limit + 1)
        if is_prime(n)
    ]


# ============================================================
# KAPPA POLYNOMIAL
# ============================================================

def f(p):

    return p * p - p + 1


def kappa(P):

    return prod(f(p) for p in P)


# ============================================================
# q-ADIC VALUATION
# ============================================================

def valuation(n, q):

    v = 0

    while n != 0 and n % q == 0:

        n //= q
        v += 1

    return v


# ============================================================
# ELEMENTARY SYMMETRIC POLYNOMIALS
# ============================================================

def elementary_symmetric(P):

    """
    Return

        e0,e1,...,ex

    where

        e0 = 1
        e1 = sum pi
        e2 = sum pi pj
        ...
        ex = product pi
    """

    e = [1]

    for p in P:

        new = [0] * (len(e) + 1)

        for j in range(len(e)):

            new[j] += e[j]
            new[j + 1] += e[j] * p

        e = new

    return e


# ============================================================
# ROOTS OF f(x) MOD m
# ============================================================

def roots_mod(m):

    return [
        x
        for x in range(m)
        if f(x) % m == 0
    ]


# ============================================================
# LOCAL LIFTING DATA
# ============================================================

def lifting_data(p, q):

    value = f(p)

    result = {

        "p": p,
        "q": q,
        "f": value,
        "valuation": valuation(value, q),

        "residues": {}
    }

    modulus = q

    for k in range(1, MAX_LIFT + 1):

        result["residues"][str(k)] = {

            "modulus": modulus,

            "p_mod":
                p % modulus,

            "f_mod":
                value % modulus,

            "root":
                value % modulus == 0,

        }

        modulus *= q

    return result


# ============================================================
# HENSEL-STYLE ROOT LIFT TEST
# ============================================================

def root_lifting_tree(q):

    """
    For each root modulo q, determine which roots survive
    modulo q^2, q^3, ...

    We don't assume Hensel's lemma blindly; we measure it.
    """

    levels = {}

    previous = roots_mod(q)

    levels["1"] = previous

    modulus = q

    for power in range(2, MAX_LIFT + 1):

        modulus *= q

        current = roots_mod(modulus)

        lifts = defaultdict(list)

        for r in current:

            parent = r % (modulus // q)

            lifts[parent].append(r)

        levels[str(power)] = {

            "roots":
                current,

            "parent_lifts":
                {
                    str(parent): children
                    for parent, children
                    in lifts.items()
                }
        }

    return levels


# ============================================================
# TEST LOCAL VALUATION DETERMINISM
# ============================================================

def valuation_residue_test(primes, q):

    """
    Check whether p mod q^k determines v_q(f(p)).

    For each k we group primes by p mod q^k and see whether
    different primes in the same residue class ever have
    different valuations.
    """

    output = {}

    for power in range(1, MAX_LIFT + 1):

        modulus = q ** power

        groups = defaultdict(set)

        for p in primes:

            groups[p % modulus].add(
                valuation(f(p), q)
            )

        ambiguous = {

            residue: sorted(values)

            for residue, values in groups.items()

            if len(values) > 1
        }

        output[str(power)] = {

            "modulus":
                modulus,

            "ambiguous_classes":
                ambiguous,

            "number_ambiguous":
                len(ambiguous),

        }

    return output


# ============================================================
# EXPAND κ AS SYMMETRIC POLYNOMIAL
# ============================================================

def symmetric_kappa_features(P):

    """
    Return elementary symmetric coordinates plus kappa.
    """

    e = elementary_symmetric(P)

    result = {

        "x": len(P),
        "primes": list(P),
        "kappa": kappa(P),
    }

    for i, value in enumerate(e):

        result[f"e{i}"] = value

    return result


# ============================================================
# SEARCH WHETHER LOWER e_i DETERMINE κ
# ============================================================

def collision_test(rows, feature_names):

    groups = defaultdict(set)

    for row in rows:

        key = tuple(
            row[name]
            for name in feature_names
        )

        groups[key].add(
            row["kappa"]
        )

    collisions = {

        key: sorted(values)

        for key, values in groups.items()

        if len(values) > 1
    }

    return collisions


# ============================================================
# MODULAR SYMMETRIC COLLISION TEST
# ============================================================

def modular_collision_test(rows, features, modulus):

    groups = defaultdict(set)

    for row in rows:

        key = tuple(
            row[name] % modulus
            for name in features
        )

        groups[key].add(
            row["kappa"] % modulus
        )

    collisions = {

        key: sorted(values)

        for key, values in groups.items()

        if len(values) > 1
    }

    return collisions


# ============================================================
# TEST κ AGAINST SIMPLE e-BASED FORMULAS
# ============================================================

def formula_tests(rows):

    results = []

    for row in rows:

        x = row["x"]

        e1 = row["e1"]
        e2 = row["e2"]

        if x >= 2:

            # Candidate suggested by the 2-body expansion.
            candidate = (
                row["e2"] ** 2
                - row["e2"] * (e1 - 1)
                + (e1 ** 2 - 2 * row["e2"])
                - e1
                + 1
            )

            results.append({

                "x":
                    x,

                "kappa":
                    row["kappa"],

                "candidate":
                    candidate,

                "error":
                    row["kappa"] - candidate,

            })

    return results


# ============================================================
# MAIN
# ============================================================

def main():

    print("=" * 72)
    print("KAPPA STRUCTURE EXPLORER III")
    print("=" * 72)

    primes = primes_up_to(MAX_PRIME)
    q_primes = [
        q
        for q in primes_up_to(Q_LIMIT)
        if q >= 3
    ]

    print(
        f"Prime limit : {MAX_PRIME}"
    )

    print(
        f"Body limit  : {MAX_BODIES}"
    )

    print(
        f"q limit     : {Q_LIMIT}"
    )

    print(
        f"Lift depth  : {MAX_LIFT}"
    )

    # ========================================================
    # A. ROOT LIFTING
    # ========================================================

    print()
    print("=" * 72)
    print("A. ROOT LIFTING")
    print("=" * 72)

    root_summary = {}

    for q in q_primes:

        roots = root_lifting_tree(q)

        root_summary[q] = roots

        print()
        print(
            f"q={q}"
        )

        print(
            "  roots mod q:",
            roots["1"]
        )

        for power in range(
            2,
            MAX_LIFT + 1
        ):

            data = roots[str(power)]

            print(
                f"  roots mod q^{power}:",
                len(data["roots"])
            )

    # ========================================================
    # B. VALUATION RESIDUE DETERMINISM
    # ========================================================

    print()
    print("=" * 72)
    print("B. VALUATION DETERMINISM")
    print("=" * 72)

    valuation_summary = {}

    for q in q_primes:

        test = valuation_residue_test(
            primes,
            q
        )

        valuation_summary[q] = test

        print()
        print(
            f"q={q}"
        )

        for power, data in test.items():

            print(
                f"  mod q^{power}: "
                f"ambiguous classes = "
                f"{data['number_ambiguous']}"
            )

    # ========================================================
    # C. GENERATE GLOBAL SYSTEMS
    # ========================================================

    print()
    print("=" * 72)
    print("C. GLOBAL SYMMETRIC STRUCTURE")
    print("=" * 72)

    rows_by_body = {}

    for x in range(
        2,
        MAX_BODIES + 1
    ):

        rows = []

        for P in combinations(
            primes,
            x
        ):

            rows.append(
                symmetric_kappa_features(P)
            )

        rows_by_body[x] = rows

        print(
            f"{x}-body systems: "
            f"{len(rows)}"
        )

    # ========================================================
    # D. EXACT SYMMETRIC COLLISIONS
    # ========================================================

    print()
    print("=" * 72)
    print("D. EXACT SYMMETRIC COLLISION TEST")
    print("=" * 72)

    exact_collision_summary = {}

    for x, rows in rows_by_body.items():

        print()
        print(
            f"{x}-body"
        )

        # progressively reveal e_i
        for number_features in range(
            1,
            x + 1
        ):

            features = [
                f"e{i}"
                for i in range(
                    1,
                    number_features + 1
                )
            ]

            collisions = collision_test(
                rows,
                features
            )

            exact_collision_summary[
                f"{x}:{features}"
            ] = len(collisions)

            print(
                f"  {features}: "
                f"{len(collisions)} collisions"
            )

    # ========================================================
    # E. MODULAR SYMMETRIC COLLISIONS
    # ========================================================

    print()
    print("=" * 72)
    print("E. MODULAR SYMMETRIC COLLISIONS")
    print("=" * 72)

    modular_summary = {}

    interesting_moduli = [
        3, 5, 7, 8,
        9, 13, 19,
        27, 49
    ]

    for x, rows in rows_by_body.items():

        modular_summary[x] = {}

        print()
        print(
            f"{x}-body"
        )

        features = [
            f"e{i}"
            for i in range(
                1,
                x + 1
            )
        ]

        for m in interesting_moduli:

            collisions = modular_collision_test(
                rows,
                features,
                m
            )

            modular_summary[x][m] = len(
                collisions
            )

            print(
                f"  mod {m:2d}: "
                f"{len(collisions)} collisions"
            )

    # ========================================================
    # F. TEST 2-BODY FORMULA
    # ========================================================

    print()
    print("=" * 72)
    print("F. TEST PAPER 2-BODY FORMULA")
    print("=" * 72)

    formula_summary = {}

    for x in range(
        2,
        MAX_BODIES + 1
    ):

        results = formula_tests(
            rows_by_body[x]
        )

        failures = [
            r
            for r in results
            if r["error"] != 0
        ]

        formula_summary[x] = {

            "tested":
                len(results),

            "failures":
                len(failures),

        }

        print(
            f"{x}-body: "
            f"tested={len(results)}, "
            f"failures={len(failures)}"
        )

    # ========================================================
    # G. SEARCH DIRECT SYMMETRIC POLYNOMIAL BEHAVIOUR
    # ========================================================

    print()
    print("=" * 72)
    print("G. κ / ELEMENTARY-SYMMETRIC RELATIONSHIPS")
    print("=" * 72)

    polynomial_observations = {}

    for x, rows in rows_by_body.items():

        observations = []

        for row in rows[:20]:

            observations.append({

                "primes":
                    row["primes"],

                "e":
                    [
                        row[f"e{i}"]
                        for i in range(x + 1)
                    ],

                "kappa":
                    row["kappa"],

            })

        polynomial_observations[x] = observations

        print()
        print(
            f"{x}-body sample:"
        )

        for observation in observations[:5]:

            print(
                observation
            )

    # ========================================================
    # H. CSV: LOCAL DATA
    # ========================================================

    print()
    print(
        f"Writing {OUTPUT_LOCAL}..."
    )

    local_fields = [
        "p",
        "q",
        "f",
        "valuation"
    ]

    for power in range(
        1,
        MAX_LIFT + 1
    ):

        local_fields += [
            f"p_mod_q{power}",
            f"f_mod_q{power}",
            f"root_mod_q{power}"
        ]

    with open(
        OUTPUT_LOCAL,
        "w",
        newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=local_fields
        )

        writer.writeheader()

        for q in q_primes:

            for p in primes:

                data = lifting_data(
                    p,
                    q
                )

                row = {

                    "p":
                        p,

                    "q":
                        q,

                    "f":
                        data["f"],

                    "valuation":
                        data["valuation"]
                }

                for power in range(
                    1,
                    MAX_LIFT + 1
                ):

                    d = data[
                        "residues"
                    ][str(power)]

                    row[
                        f"p_mod_q{power}"
                    ] = d["p_mod"]

                    row[
                        f"f_mod_q{power}"
                    ] = d["f_mod"]

                    row[
                        f"root_mod_q{power}"
                    ] = d["root"]

                writer.writerow(row)

    # ========================================================
    # I. CSV: GLOBAL DATA
    # ========================================================

    print(
        f"Writing {OUTPUT_GLOBAL}..."
    )

    global_fields = [
        "x",
        "primes",
        "kappa"
    ]

    for i in range(
        0,
        MAX_BODIES + 1
    ):

        global_fields.append(
            f"e{i}"
        )

    with open(
        OUTPUT_GLOBAL,
        "w",
        newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=global_fields
        )

        writer.writeheader()

        for x, rows in rows_by_body.items():

            for row in rows:

                output = {

                    "x":
                        row["x"],

                    "primes":
                        str(row["primes"]),

                    "kappa":
                        row["kappa"],
                }

                for i in range(
                    0,
                    MAX_BODIES + 1
                ):

                    output[
                        f"e{i}"
                    ] = row.get(
                        f"e{i}",
                        ""
                    )

                writer.writerow(output)

    # ========================================================
    # J. SUMMARY JSON
    # ========================================================

    summary = {

        "configuration": {

            "MAX_PRIME":
                MAX_PRIME,

            "MAX_BODIES":
                MAX_BODIES,

            "Q_LIMIT":
                Q_LIMIT,

            "MAX_LIFT":
                MAX_LIFT,

        },

        "root_lifting":
            root_summary,

        "valuation_determinism":
            valuation_summary,

        "exact_symmetric_collisions":
            exact_collision_summary,

        "modular_symmetric_collisions":
            modular_summary,

        "paper_formula":
            formula_summary,

        "polynomial_observations":
            polynomial_observations,

    }

    print(
        f"Writing {OUTPUT_SUMMARY}..."
    )

    with open(
        OUTPUT_SUMMARY,
        "w"
    ) as file:

        json.dump(
            summary,
            file,
            indent=2
        )

    # ========================================================
    # DONE
    # ========================================================

    print()
    print("=" * 72)
    print("EXPLORATION III COMPLETE")
    print("=" * 72)

    print()
    print("Send me:")
    print(
        f"  {OUTPUT_SUMMARY}"
    )
    print(
        f"  {OUTPUT_LOCAL}"
    )
    print(
        f"  {OUTPUT_GLOBAL}"
    )

    print()
    print(
        "JSON first. The two CSV files are useful if "
        "we find an interesting structure."
    )


if __name__ == "__main__":
    main()