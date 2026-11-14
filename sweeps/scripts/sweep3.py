"""Sweep 3: which L-values / products have BBP formulas in base p^m, length n?
Usage: python -I sweep3.py p s "m1:n1,m2:n2,..." "target1;target2;..." [label]
Targets use targets2.values syntax (e.g. s2*L(-8,2), pi*l2)."""
import sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import mpmath as mp
import qvf, targets2

def run(p, s, m, n, names):
    d = len(names) + n
    dps = max(600, 45 * d)
    t0 = time.time()
    b = p ** m
    T = qvf.bbp_basis(s, b, n, dps)
    tv = targets2.values(names, dps)
    vals = tv + T
    digits = dps - 20
    rels, dt = qvf.lll_relations(vals, digits)
    goods = qvf.good_relations(rels, digits)
    k = len(names)
    pure = [g for g in goods if not any(g[k:])]   # relations among targets only (must be none)
    zero = [g for g in goods if not any(g[:k])]
    mixed = [g for g in goods if any(g[:k]) and any(g[k:])]
    rank, rows = qvf.proj_rank(mixed + pure, k)
    print(f"=== p={p} s={s} b={p}^{m} n={n} dim={d} dps={dps} | rel {len(goods)} (zero {len(zero)}, target-only {len(pure)}) | target rank {rank} | {time.time()-t0:.1f}s", flush=True)
    for r in rows:
        print("   in span:", ", ".join(f"{nm}:{x}" for nm, x in zip(names, r) if x != 0), flush=True)
    for g in pure:
        print("   WARNING target-only relation:", g[:k], flush=True)

def main():
    p, s = int(sys.argv[1]), int(sys.argv[2])
    cfgs = [tuple(int(x) for x in c.split(':')) for c in sys.argv[3].split(',')]
    names = sys.argv[4].split(';')
    for m, n in cfgs:
        run(p, s, m, n, names)
    print("SWEEP3 DONE", flush=True)

if __name__ == '__main__':
    main()
