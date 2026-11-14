"""Binary BBP span search with conductor-5/20/40 targets.  Usage: python -I span40.py m n [digits_per_dim]"""
import sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0] + '/../../2026-11-11_L8_binary_BBP/scripts')
import mpmath as mp
import qvf, targets2

def tvals(dps):
    mp.mp.dps = dps + 30
    pi, l2 = mp.pi, mp.log(2)
    s2, s5, s10 = mp.sqrt(2), mp.sqrt(5), mp.sqrt(10)
    G = mp.catalan
    L8 = targets2.Lval(-8, 2); L20 = targets2.Lval(-20, 2); L40 = targets2.Lval(-40, 2)
    mp.mp.dps = dps + 30
    C1 = mp.clsin(2, 2 * pi / 5); C2 = mp.clsin(2, 4 * pi / 5)
    a5 = mp.sqrt(10 + 2 * s5); b5 = mp.sqrt(10 - 2 * s5)
    T = [('G', G), ('L20', L20), ('s5*L20', s5 * L20),
         ('C1', C1), ('C2', C2), ('s5*C1', s5 * C1), ('s5*C2', s5 * C2), ('a5*C1', a5 * C1), ('a5*C2', a5 * C2), ('b5*C1', b5 * C1), ('b5*C2', b5 * C2),
         ('pi*l2', pi * l2), ('s5*pi*l2', s5 * pi * l2), ('a5*pi*l2', a5 * pi * l2), ('b5*pi*l2', b5 * pi * l2),
         ('pi^2', pi ** 2), ('l2^2', l2 ** 2)]
    return [t[0] for t in T], [t[1] for t in T]

def main():
    m, n = int(sys.argv[1]), int(sys.argv[2])
    per = int(sys.argv[3]) if len(sys.argv) > 3 else 45
    names, _ = tvals(50)
    d = len(names) + n
    dps = max(600, per * d)
    t0 = time.time()
    T = qvf.bbp_basis(2, 2 ** m, n, dps)
    names, tv = tvals(dps)
    vals = tv + T
    digits = dps - 20
    rels, dt = qvf.lll_relations(vals, digits)
    goods = qvf.good_relations(rels, digits)
    k = len(names)
    rank, rows = qvf.proj_rank(goods, k)
    print(f"=== b=2^{m} n={n} dim={d} dps={dps} | rel {len(goods)} | target rank {rank} | LLL {dt:.1f}s total {time.time()-t0:.1f}s", flush=True)
    for r in rows:
        print("   in span:", ", ".join(f"{nm}:{x}" for nm, x in zip(names, r) if x != 0), flush=True)
    pure = [g for g in goods if not any(g[k:])]
    for g in pure:
        print("   WARNING target-only relation:", g[:k])

if __name__ == '__main__':
    main()
