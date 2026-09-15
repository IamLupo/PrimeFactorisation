#!/usr/bin/env python3
import csv, itertools, math
from collections import defaultdict
from pathlib import Path

OUT = Path("kappa_aux_results")
OUT.mkdir(exist_ok=True)

def F(x):
    return x*x - x + 1

def primes_upto(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p*p:limit+1:p] = b"\x00" * (((limit-p*p)//p) + 1)
    return [i for i in range(2, limit + 1) if sieve[i]]

def bits(x):
    return x.bit_length()

def A_of(p, q):
    return F(p) * F(q)

def B_of(rs):
    B = 1
    for r in rs:
        B *= F(r)
    return B

def invariants(p, q):
    n = p*q
    u = p+q
    A = A_of(p,q)
    D = 4*A - 3*(n-1)**2
    return n,u,A,D

def K_of(A, rs):
    return 1 - A*B_of(rs)

def recover_A(K, rs):
    B = B_of(rs)
    z = 1-K
    return z//B if z % B == 0 else None

def recover_factors(n, A):
    D = 4*A - 3*(n-1)**2
    if D < 0:
        return []
    t = math.isqrt(D)
    if t*t != D:
        return []
    out = []
    for sign in (-1, 1):
        z = n+1+sign*t
        if z % 2:
            continue
        u = z//2
        d = u*u-4*n
        if d < 0:
            continue
        sd = math.isqrt(d)
        if sd*sd == d:
            a,b = u-sd, u+sd
            if a*b == n:
                out.append(tuple(sorted((a,b))))
    return sorted(set(out))

def sequences(primes, max_terms=15, target_bits=100):
    seq = []
    candidates = [
        ("first", primes[:max_terms]),
        ("odd_first", [p for p in primes if p != 2][:max_terms]),
        ("every_other", primes[::2][:max_terms]),
    ]
    B = 1
    rs = []
    for r in primes:
        rs.append(r)
        B *= F(r)
        if bits(B) >= target_bits:
            break
    candidates.append(("to_target_bits", rs))
    seen = set()
    for name, xs in candidates:
        t = tuple(xs)
        if t not in seen:
            seen.add(t)
            seq.append((name,t))
    return seq

def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

def main():
    primes = primes_upto(500)
    pairs = list(itertools.combinations(primes,2))[:2500]
    seqs = sequences(primes)

    # 1. How fast does B=product F(r) reach 100 bits?
    rows = []
    B = 1
    for i,r in enumerate(primes[:50],1):
        B *= F(r)
        rows.append({
            "i":i, "r":r, "F_r":F(r), "B":B,
            "B_bits":bits(B), "B_digits":len(str(B))
        })
    write_csv(OUT/"auxiliary_growth.csv", rows)

    # 2. Exact controlled experiments.
    rows = []
    for p,q in pairs:
        n,u,A,D = invariants(p,q)
        for name,rs in seqs:
            B = B_of(rs)
            K = 1-A*B
            A2 = recover_A(K,rs)
            cand = recover_factors(n,A2) if A2 is not None else []
            rows.append({
                "p":p,"q":q,"n":n,"n_bits":bits(n),"A":A,
                "aux":name,"rs":",".join(map(str,rs)),
                "r_count":len(rs),"B":B,"B_bits":bits(B),
                "K":K,"recovered_A":A2,
                "A_recovered":A2==A,
                "factor_recovered":(p,q) in cand,
                "n_times_B_bits":bits(n*B)
            })
    write_csv(OUT/"exact_auxiliary_experiments.csv",rows)

    # 3. Test whether modular observations become more informative.
    moduli = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
    rows = []
    for name,rs in seqs:
        B = B_of(rs)
        for m in moduli:
            buckets = defaultdict(set)
            for p,q in pairs:
                n,u,A,D = invariants(p,q)
                K = 1-A*B
                # Public quantities n,m,B; observed quantity K mod m.
                buckets[(n%m,K%m)].add((p,q))
            rows.append({
                "aux":name,"r_count":len(rs),"B_bits":bits(B),
                "modulus":m,
                "signatures":len(buckets),
                "unique":sum(len(v)==1 for v in buckets.values()),
                "max_collision":max(map(len,buckets.values()))
            })
    write_csv(OUT/"modular_information.csv",rows)

    # 4. GCD/divisibility behavior of B against n and simple n-polynomials.
    rows = []
    for p,q in pairs[:1000]:
        n,u,A,D = invariants(p,q)
        for name,rs in seqs:
            B = B_of(rs)
            rows.append({
                "p":p,"q":q,"n":n,"aux":name,"r_count":len(rs),
                "B_bits":bits(B),
                "B_mod_n":B%n,
                "B_mod_n_minus_1":B%(n-1),
                "B_mod_n_plus_1":B%(n+1),
                "gcd_B_n":math.gcd(B,n),
                "gcd_B_n_minus_1":math.gcd(B,n-1),
                "gcd_B_n_plus_1":math.gcd(B,n+1),
                "gcd_B_n2_minus_1":math.gcd(B,n*n-1),
                "gcd_A_B":math.gcd(A,B)
            })
    write_csv(OUT/"B_gcd_structure.csv",rows)

    # 5. For fixed n, compare collisions of A, K and K modulo m.
    #    This is the most direct "does adding r's narrow the space?"
    by_n = defaultdict(list)
    for p,q in pairs:
        n,u,A,D = invariants(p,q)
        by_n[n].append((p,q,A))

    rows = []
    for name,rs in seqs:
        B = B_of(rs)
        for n,items in by_n.items():
            if len(items) < 2:
                continue
            exact_A = defaultdict(list)
            exact_K = defaultdict(list)
            for p,q,A in items:
                exact_A[A].append((p,q))
                exact_K[1-A*B].append((p,q))
            for m in moduli:
                red = defaultdict(list)
                for p,q,A in items:
                    red[(1-A*B)%m].append((p,q))
                rows.append({
                    "n":n,"aux":name,"r_count":len(rs),
                    "B_bits":bits(B),"factor_pairs":len(items),
                    "A_values":len(exact_A),
                    "K_values":len(exact_K),
                    "K_mod_m_values":len(red),
                    "modulus":m,
                    "K_mod_collisions":sum(len(v)>1 for v in red.values())
                })
    write_csv(OUT/"fixed_n_collision_analysis.csv",rows)

    print("="*78)
    print("KAPPA AUXILIARY-PRIME / INFORMATION-GROWTH EXPERIMENT")
    print("="*78)
    print()
    print("Identity:")
    print("  K = 1 - A*B")
    print("  A = F(p)F(q)")
    print("  B = product F(r_i)")
    print()
    print("If r_i and K are known, A=(1-K)/B exactly.")
    print("The experiment therefore asks whether many auxiliary cases")
    print("produce extra constraints that can be obtained from n alone.")
    print()
    print("Auxiliary sequences:")
    for name,rs in seqs:
        B=B_of(rs)
        print(f"  {name:14s}: {len(rs):2d} primes, B bits={bits(B):3d}, rs={rs}")
    print()
    print("Example hidden targets:")
    for p,q in pairs[:8]:
        n,u,A,D=invariants(p,q)
        print(f"  p={p:3d} q={q:3d} n={n:6d} A={A:8d} n_bits={bits(n)}")
    print()
    print("Files:")
    for x in [
        "auxiliary_growth.csv",
        "exact_auxiliary_experiments.csv",
        "modular_information.csv",
        "B_gcd_structure.csv",
        "fixed_n_collision_analysis.csv",
    ]:
        print(" ",OUT/x)
    print()
    print("Interpretation to look for:")
    print("  1. B reaches 100 bits using surprisingly few small primes.")
    print("  2. Exact K always recovers A; auxiliary factors do not change that.")
    print("  3. The interesting result would be a new n-only restriction on A.")
    print("  4. Pay particular attention to fixed_n_collision_analysis.csv.")
    print("     If K residues reduce collisions for many moduli, record that")
    print("     as evidence of information, but not as a factoring algorithm.")
    print()
    print("DONE")

if __name__ == "__main__":
    main()
