"""Extract explicit binary BBP formulas C = (1/c0) P(2, 2^m, n, A) for the conductor-5 / conductor-20 constants
and verify at a higher precision.  Usage: python -I extract40.py which m n verify_dps
which: L20 | S1 | S2 (constants defined below)."""
import sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/../../2026-11-11_L8_binary_BBP/scripts')
import mpmath as mp
import qvf
from math import gcd

def hurwitz_L(chi, q, s=2):
    # sum_{n>=1} chi(n)/n^s with chi periodic mod q
    return mp.fsum(chi[r % q] * mp.zeta(s, mp.mpf(r) / q) for r in range(1, q + 1) if chi[r % q]) / mp.mpf(q) ** s

def const(which, dps):
    mp.mp.dps = dps + 40
    pi, l2, s5 = mp.pi, mp.log(2), mp.sqrt(5)
    al, be = mp.sqrt(10 + 2 * s5), mp.sqrt(10 - 2 * s5)
    chi20 = {r: 0 for r in range(20)}
    for r in (1, 3, 7, 9): chi20[r] = 1
    for r in (11, 13, 17, 19): chi20[r] = -1
    if which == 'L20':
        return hurwitz_L(chi20, 20) - s5 / 15 * pi * l2
    S1 = hurwitz_L({0: 0, 1: 1, 2: 0, 3: 0, 4: -1}, 5)
    S2 = hurwitz_L({0: 0, 1: 0, 2: 1, 3: -1, 4: 0}, 5)
    if which == 'S1':   # 5 S1 - (pi log2/1105)(128 al + 36 be)
        return 5 * S1 - pi * l2 / 1105 * (128 * al + 36 * be)
    if which == 'S2':   # 5 S2 + (pi log2/1105)(36 al - 128 be)
        return 5 * S2 + pi * l2 / 1105 * (36 * al - 128 * be)
    raise ValueError

def main():
    which = sys.argv[1]; m, n = int(sys.argv[2]), int(sys.argv[3]); vdps = int(sys.argv[4])
    b = 2 ** m
    d = n + 1
    dps = max(800, 50 * d)
    T = qvf.bbp_basis(2, b, n, dps)
    C = const(which, dps)
    rels, dt = qvf.lll_relations([C] + T, dps - 20)
    goods = qvf.good_relations(rels, dps - 20)
    withC = [g for g in goods if g[0] != 0]
    print(f"{which}: b=2^{m} n={n} dps={dps}: relations {len(goods)}, with C {len(withC)}, LLL {dt:.1f}s", flush=True)
    if not withC:
        print("NO FORMULA"); return
    withC.sort(key=lambda g: sum(abs(x) for x in g))
    g = withC[0]
    c0 = g[0]; A = [-x for x in g[1:]]
    if c0 < 0: c0, A = -c0, [-a for a in A]
    gg = 0
    for a in A + [c0]: gg = gcd(gg, a)
    c0 //= gg; A = [a // gg for a in A]
    print(f"FORMULA {which} = (1/{c0}) * P(2, 2^{m}, {n}, {A})", flush=True)
    t0 = time.time()
    Tv = qvf.bbp_basis(2, b, n, vdps)
    Cv = const(which, vdps)
    diff = Cv - mp.fsum(mp.mpf(a) * t for a, t in zip(A, Tv)) / c0
    print(f"verify at {vdps} digits: log10|LHS-RHS| = {float(mp.log10(abs(diff))) if diff != 0 else '-inf'} ({time.time()-t0:.1f}s)", flush=True)

if __name__ == '__main__':
    main()
