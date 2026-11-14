"""Extract an explicit BBP formula C = (1/r) P(s, b, n, A) for a single constant C, then verify it
at a higher precision.  Usage: python -I extract.py p s m n "<target expr>" [verify_dps]
Target syntax as in targets2.values (one constant)."""
import sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import mpmath as mp
import qvf, targets2
from fractions import Fraction
from math import gcd

def combo(expr, dps):
    # expr: 'q1:term1;q2:term2;...' (rational q, target syntax of targets2) or a single term
    if ':' not in expr:
        return targets2.values([expr], dps)[0]
    tot = mp.mpf(0)
    for part in expr.split(';'):
        q, t = part.split(':', 1)
        tot += mp.mpf(Fraction(q).numerator) / Fraction(q).denominator * targets2.values([t], dps)[0]
    return tot

def main():
    p, s, m, n = (int(x) for x in sys.argv[1:5])
    expr = sys.argv[5]
    vdps = int(sys.argv[6]) if len(sys.argv) > 6 else 4000
    b = p ** m
    d = n + 1
    dps = max(800, 50 * d)
    T = qvf.bbp_basis(s, b, n, dps)
    C = combo(expr, dps)
    rels, dt = qvf.lll_relations([C] + T, dps - 20)
    goods = qvf.good_relations(rels, dps - 20)
    withC = [g for g in goods if g[0] != 0]
    zero = [g for g in goods if g[0] == 0]
    print(f"config p={p} s={s} b={p}^{m} n={n}: relations {len(goods)} (zero relations {len(zero)}) LLL {dt:.1f}s")
    if not withC:
        print("NO FORMULA"); return
    # smallest relation involving C
    withC.sort(key=lambda g: sum(abs(x) for x in g))
    g = withC[0]
    c0 = g[0]
    A = [-x for x in g[1:]]
    if c0 < 0:
        c0, A = -c0, [-a for a in A]
    gg = 0
    for a in A + [c0]:
        gg = gcd(gg, a)
    c0 //= gg; A = [a // gg for a in A]
    print(f"FORMULA: {expr} = (1/{c0}) * P({s}, {p}^{m}, {n}, {A})")
    # verify at higher precision
    t0 = time.time()
    Tv = qvf.bbp_basis(s, b, n, vdps)
    Cv = combo(expr, vdps)
    rhs = mp.fsum(mp.mpf(a) * t for a, t in zip(A, Tv)) / c0
    diff = Cv - rhs
    print(f"verify at {vdps} digits: |LHS-RHS| = {mp.nstr(abs(diff), 5)} (log10 {float(mp.log10(abs(diff))) if diff != 0 else '-inf'}), {time.time()-t0:.1f}s")
    for z in zero[:3]:
        print("  zero relation:", z[1:])

if __name__ == '__main__':
    main()
