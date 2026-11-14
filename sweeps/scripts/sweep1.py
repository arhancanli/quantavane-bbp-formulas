"""Sweep 1: ternary (base 3^m) BBP spans, degrees 2 and 3, many lengths.

Question: which target constants (or combinations) lie in the Q-span of
P(s, 3^m, n, A)?  Main aim: a third independent relation among
{zeta(3), pi^2 log 3, log^3 3}, which would give a ternary formula for zeta(3) alone.
Writes one block per configuration to stdout (run detached with a log).
"""
import sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import mpmath as mp
import qvf
from exp1_span import target_values

T3 = ['z3', 'pi2l3', 'l3_3', 'pi2l2', 'l2_3', 'l2l3_2', 'l2_2l3',
      's3pi3', 's3pil3_2', 's3pil2_2', 's3pil2l3', 's3L12_3', 'L12_3',
      'pi3', 'pil3_2', 'pil2_2', 'pil2l3', 'Gl2', 'Gl3']
T2 = ['pi2', 'l3_2', 's3pil3', 'pil3', 'l2_2', 'l2l3', 'pil2', 's3pil2',
      's3Lm3_2', 'G']


def run(s, b, n, names):
    d = len(names) + n
    dps = max(600, 45 * d)
    t0 = time.time()
    T = qvf.bbp_basis(s, b, n, dps)
    tv = target_values(names, dps)
    vals = tv + T
    digits = dps - 20
    rels, dt = qvf.lll_relations(vals, digits)
    goods = qvf.good_relations(rels, digits)
    k = len(names)
    rank, rows = qvf.proj_rank(goods, k)
    zero_rel = sum(1 for g in goods if not any(g[:k]))
    print(f"=== s={s} b=3^{round(mp.log(b,3))} n={n} dim={d} dps={dps} | relations {len(goods)} (zero relations {zero_rel}) | target rank {rank} | {time.time()-t0:.1f}s", flush=True)
    for r in rows:
        print("   in span:", ", ".join(f"{nm}:{x}" for nm, x in zip(names, r) if x != 0), flush=True)


def main():
    configs = []
    for n in [4, 6, 8, 10, 12, 16, 18, 20, 24]:
        for frac_num, frac_den in [(1, 4), (1, 3), (1, 2), (2, 3), (1, 1)]:
            if (n * frac_num) % frac_den:
                continue
            m = n * frac_num // frac_den
            if m < 1:
                continue
            configs.append((m, n))
    configs = sorted(set(configs), key=lambda c: (c[1], c[0]))
    which = sys.argv[1] if len(sys.argv) > 1 else '3'
    for (m, n) in configs:
        b = 3 ** m
        if which in ('2', 'both'):
            run(2, b, n, T2)
        if which in ('3', 'both'):
            run(3, b, n, T3)
    print("SWEEP1 DONE", flush=True)


if __name__ == '__main__':
    main()
