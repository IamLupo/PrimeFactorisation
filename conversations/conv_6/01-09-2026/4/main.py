#!/usr/bin/env python3
from __future__ import annotations
import math, random, time
from collections import defaultdict
from typing import Dict, List, Sequence, Tuple

SEED=1511464998
SCALES=(10**9,10**12,10**16)
ANCHORS_PER_SCALE=2
TARGET_T=(1111,1000,909)
R_CASES_PER_N=3
AUX_COUNT=4
AUX_X_MAX=500_000
STATIC_K_MAX=5000
AUX_K_WINDOW=24
SMALL_PRIMES=(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997)

def is_prime(n:int)->bool:
    if n<2:return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%p==0:return n==p
    d=n-1;s=0
    while d%2==0:s+=1;d//=2
    for a in (2,325,9375,28178,450775,9780504,1795265022):
        if a%n==0:continue
        x=pow(a,d,n)
        if x in (1,n-1):continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1:break
        else:return False
    return True

def next_prime(n:int)->int:
    x=max(2,n)
    if x==2:return 2
    if x%2==0:x+=1
    while not is_prime(x):x+=2
    return x

def prev_prime(n:int)->int:
    if n<=2:return 2
    x=n if n%2 else n-1
    while x>2 and not is_prime(x):x-=2
    return x

def build_static_k_table(k_max:int):
    table=defaultdict(list)
    for K in range(1,k_max+1):
        for d in range(1,math.isqrt(K)+1):
            if K%d==0:
                e=K//d; table[K].append((d,e))
                if d!=e: table[K].append((e,d))
    return {k:tuple(v) for k,v in table.items()}

def generate_semiprime(scale:int,rng:random.Random):
    root=math.isqrt(scale)
    while True:
        p=next_prime(rng.randint(max(1000,int(root*.72)),int(root*1.35)))
        q=next_prime(rng.randint(max(1000,int(root*.72)),int(root*1.35)))
        if p!=q:return (p*q,max(p,q),min(p,q))

def choose_moduli(n:int,target_t:int,case_index:int):
    base=math.isqrt(max(3,n//target_t))
    tweaks=((-23,19),(-7,11),(17,-13))
    d1,d2=tweaks[case_index%len(tweaks)]
    r1=prev_prime(base+d1);r2=next_prime(base+d2)
    if r1==r2:r2=next_prime(r2+2)
    return r1,r2

def coords(n,r1,r2,p,q):
    k,a=divmod(p,r1);l,b=divmod(q,r2);T=n//(r1*r2);K=k*l
    return dict(k=k,l=l,a=a,b=b,T=T,K=K,E=T-K)

def find_easy_neighbors(n:int,count:int):
    found=[]
    for x in range(1,AUX_X_MAX+1,2):
        N=n+x
        for d in SMALL_PRIMES:
            if N%d==0 and is_prime(N//d):
                found.append(dict(x=x,N=N,px=d,qx=N//d))
                break
        if len(found)>=count:break
    if len(found)<count:raise RuntimeError(f'Could not find {count} easy n+x neighbors up to x={AUX_X_MAX}')
    return found

def aux_records(neighbors,r1,r2):
    out=[]
    for z in neighbors:
        for px,qx,ori in ((z['px'],z['qx'],0),(z['qx'],z['px'],1)):
            kx,ax=divmod(px,r1);lx,bx=divmod(qx,r2);Tx=z['N']//(r1*r2);Kx=kx*lx;Ex=Tx-Kx
            out.append(dict(x=z['x'],N=z['N'],px=px,qx=qx,orientation=ori,kx=kx,lx=lx,ax=ax,bx=bx,Tx=Tx,Kx=Kx,Ex=Ex))
    return out

def carry_interval(c,mult,mod):
    if mult<=0:return (1,0)
    lo=(c*mod+mult-1)//mult
    hi=(((c+1)*mod+mult-1)//mult)-1
    return lo,hi

def c3_exact(a,b,k,l,r1,r2):
    c1=(a*l)//r1;c2=(b*k)//r2
    d1=a*l-c1*r1;d2=b*k-c2*r2
    return (d1*r2+d2*r1+a*b)//(r1*r2)

def reconstruct_ab(n,r1,r2,k,l,E):
    if k<=0 or l<=0:return []
    A=k*r1;B=l*r2;hits=[]
    for c3 in (0,1,2):
        rem=E-c3
        if rem<0:continue
        c1_lo=max(0,rem-(k-1));c1_hi=min(l-1,rem)
        for c1 in range(c1_lo,c1_hi+1):
            c2=rem-c1
            if not 0<=c2<k:continue
            al,ah=carry_interval(c1,l,r1);bl,bh=carry_interval(c2,k,r2)
            al=max(al,0);ah=min(ah,r1-1);bl=max(bl,0);bh=min(bh,r2-1)
            if al>ah or bl>bh:continue
            ah=min(ah,n//(B+bl)-A)
            al=max(al,(n+(B+bh)-1)//(B+bh)-A)
            if al>ah:continue
            for a in range(al,ah+1):
                den=A+a
                if den<=0 or n%den:continue
                b=n//den-B
                if not (bl<=b<=bh and 0<=b<r2):continue
                if c3_exact(a,b,k,l,r1,r2)!=c3:continue
                if den*(B+b)!=n:continue
                hits.append(dict(k=k,l=l,a=a,b=b,c1=c1,c2=c2,c3=c3,carry_sum=c1+c2+c3,a_width=ah-al+1,b_width=bh-bl+1))
                if len(hits)>=8:return hits
    return hits

def aux_seed_k_values(records,kmax):
    s=set()
    for r in records:
        for d in range(-AUX_K_WINDOW,AUX_K_WINDOW+1):
            K=r['Kx']+d
            if 1<=K<=kmax:s.add(K)
    return sorted(s)

def reconstruct_from_k(n,r1,r2,T,Kcands,table):
    out=[];seen=set()
    for K in Kcands:
        E=T-K
        if E<0:continue
        for k,l in table.get(K,()):
            for h in reconstruct_ab(n,r1,r2,k,l,E):
                key=(K,k,l,h['a'],h['b'])
                if key not in seen:
                    seen.add(key);o=dict(h);o.update(K=K,E=E);out.append(o)
    return out

def true_present(sols,t):
    return any(s['K']==t['K'] and s['k']==t['k'] and s['l']==t['l'] and s['a']==t['a'] and s['b']==t['b'] for s in sols)

def fmt(x):return f'{x:,}'

def main():
    rng=random.Random(SEED);t0=time.perf_counter();table=build_static_k_table(STATIC_K_MAX);build=time.perf_counter()-t0
    print('='*100);print('START EXPERIMENT 78');print('AUXILIARY n+x -> Kx/Ex/Tx/(kx,lx)/(ax,bx) -> K -> (k,l) -> (a,b)');print('='*100)
    print(f'configuration: scales={[f"1e{len(str(s))-1}" for s in SCALES]}, anchors={ANCHORS_PER_SCALE}, Rcases={R_CASES_PER_N}, aux={AUX_COUNT}, Kmax={STATIC_K_MAX}, auxKwindow=+/-{AUX_K_WINDOW}, seed={SEED}')
    total=aux_total=exact_hits=seed_hits=aux_rec=fall_rec=idfails=0
    for scale in SCALES:
        print('\n'+'='*100);print(f'SCALE 1e+{len(str(scale))-1}');print('='*100)
        for ai in range(ANCHORS_PER_SCALE):
            n,p,q=generate_semiprime(scale,rng);print(f'\nANCHOR {ai+1}/{ANCHORS_PER_SCALE} n={fmt(n)}');print(f'    generator factors = ({fmt(p)},{fmt(q)}) [validation only]')
            for ci in range(R_CASES_PER_N):
                total+=1;r1,r2=choose_moduli(n,TARGET_T[ci],ci);t=coords(n,r1,r2,p,q);neigh=find_easy_neighbors(n,AUX_COUNT);recs=aux_records(neigh,r1,r2);aux_total+=len(recs)
                ek=sum(r['Kx']==t['K'] for r in recs);exact_hits+=ek;seeds=aux_seed_k_values(recs,STATIC_K_MAX);inside=t['K'] in seeds;seed_hits+=inside
                print(f'\nCASE {ci+1}/{R_CASES_PER_N}: r1={fmt(r1)} r2={fmt(r2)} T={t["T"]} TRUE K={t["K"]} E={t["E"]} (k,l)=({t["k"]},{t["l"]}) (a,b)=({fmt(t["a"])},{fmt(t["b"])})')
                print(f'    aux records={len(recs)} exact Kx==K={ek} aux-seed-K count={len(seeds)} true-K-in-seeds={inside}')
                for r in recs:
                    print(f'        x={fmt(r["x"]):>7} N={fmt(r["N"])} (px,qx)=({fmt(r["px"])},{fmt(r["qx"])}) Tx={r["Tx"]:>4} Ex={r["Ex"]:>4} Kx={r["Kx"]:>4} (kx,lx)=({fmt(r["kx"])},{fmt(r["lx"])}) (ax,bx)=({fmt(r["ax"])},{fmt(r["bx"])})')
                fails=0
                for r in recs:
                    if r['Ex']-(t['E']) != (r['Tx']-t['T'])-(r['Kx']-t['K']):fails+=1
                idfails+=fails;print(f'    exact auxiliary identity failures={fails}/{len(recs)}')
                sols=reconstruct_from_k(n,r1,r2,t['T'],seeds,table);hit=true_present(sols,t);aux_rec+=hit
                print(f'    AUX-SEED RECOVERY: solutions={len(sols)} true-recovered={hit}')
                for s in sols[:5]:print(f'        K={s["K"]} E={s["E"]} (k,l)=({s["k"]},{s["l"]}) (a,b)=({fmt(s["a"])},{fmt(s["b"])}) c=({s["c1"]},{s["c2"]},{s["c3"]}) a-search-width={fmt(s["a_width"])}')
                if not hit:
                    sols2=reconstruct_from_k(n,r1,r2,t['T'],range(1,STATIC_K_MAX+1),table);hit2=true_present(sols2,t);fall_rec+=hit2
                    print(f'    STATIC FALLBACK: solutions={len(sols2)} true-recovered={hit2}')
                    for s in sols2[:5]:print(f'        K={s["K"]} E={s["E"]} (k,l)=({s["k"]},{s["l"]}) (a,b)=({fmt(s["a"])},{fmt(s["b"])} ) c=({s["c1"]},{s["c2"]},{s["c3"]}) a-search-width={fmt(s["a_width"])}')
                c1=(t['a']*t['l'])//r1;c2=(t['b']*t['k'])//r2;c3=c3_exact(t['a'],t['b'],t['k'],t['l'],r1,r2)
                print(f'    TRUE CARRY STATE: c=({c1},{c2},{c3}) sum={c1+c2+c3} E={t["E"]}')
    elapsed=time.perf_counter()-t0
    print('\n'+'='*100);print('GLOBAL SUMMARY');print('='*100)
    print(f'    build time={build:.6f}s');print(f'    cases={total} auxiliary records={aux_total}');print(f'    exact Kx==trueK records={exact_hits}');print(f'    true K in aux-seed sets={seed_hits}/{total}');print(f'    aux-seed true recoveries={aux_rec}/{total}');print(f'    static fallback true recoveries={fall_rec}/{total}');print(f'    auxiliary identity failures={idfails}');print(f'    total runtime={elapsed:.4f}s')
    print('\n'+'='*100);print('KEY QUESTION');print('='*100)
    print('    The reconstruction stage is: K -> (k,l) -> E=T-K -> carry states -> (a,b).')
    print('    Auxiliary observations provide Kx, Ex, Tx and residue coordinates.')
    print('    The exact identity Ex-E=(Tx-T)-(Kx-K) is checked explicitly.')
    print('    AUX-SEED RECOVERY means K was found from the auxiliary Kx neighborhood without scanning the full static K table.')
    print('    STATIC FALLBACK RECOVERY means the carry reconstruction succeeds once K is allowed to come from the complete static table.')
    print('    No hidden p,q are used during reconstruction; they are used only for final validation.')
    print('='*100);print('FINISHED EXPERIMENT 78');print('='*100)

if __name__=='__main__':main()
