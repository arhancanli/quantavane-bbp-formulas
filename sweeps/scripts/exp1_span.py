"""Experiment 1: which target constants lie in the Q-span of base-b, length-n, degree-s BBP sums?

Usage: python -I exp1_span.py s b n dps target1,target2,...
Prints the LLL relation lattice projected onto the targets.
"""
import sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import mpmath as mp
import qvf


def target_values(names, dps):
    mp.mp.dps = dps + 30
    pi, l2, l3 = mp.pi, mp.log(2), mp.log(3)
    s3 = mp.sqrt(3)
    z3, z5 = mp.zeta(3), mp.zeta(5)
    G = mp.catalan
    def L(dchar, s):
        # Dirichlet L-value for a real primitive character given as a list over residues mod m
        m = len(dchar)
        return mp.fsum(dchar[r] * mp.zeta(s, mp.mpf(r) / m) for r in range(1, m) if dchar[r]) / mp.mpf(m) ** s
    chi12 = [0, 1, 0, 0, 0, -1, 0, -1, 0, 0, 0, 1]
    chim3 = [0, 1, -1]
    chim4 = [0, 1, 0, -1]
    chi5 = [0, 1, -1, -1, 1]
    chim8 = [0, 1, 0, 1, 0, -1, 0, -1]
    chi8 = [0, 1, 0, -1, 0, -1, 0, 1]
    table = {
        'one': lambda: mp.mpf(1),
        'z3': lambda: z3, 'z5': lambda: z5,
        'pi': lambda: pi, 'pi2': lambda: pi ** 2, 'pi3': lambda: pi ** 3, 'pi4': lambda: pi ** 4,
        's3pi': lambda: s3 * pi, 's3pi3': lambda: s3 * pi ** 3,
        'l2': lambda: l2, 'l3': lambda: l3, 'l5': lambda: mp.log(5), 'l7': lambda: mp.log(7), 'l13': lambda: mp.log(13),
        'l2_2': lambda: l2 ** 2, 'l3_2': lambda: l3 ** 2, 'l2l3': lambda: l2 * l3,
        'pil2': lambda: pi * l2, 'pil3': lambda: pi * l3, 's3pil3': lambda: s3 * pi * l3, 's3pil2': lambda: s3 * pi * l2,
        'pi2l3': lambda: pi ** 2 * l3, 'pi2l2': lambda: pi ** 2 * l2,
        'l3_3': lambda: l3 ** 3, 'l2_3': lambda: l2 ** 3, 'l2l3_2': lambda: l2 * l3 ** 2, 'l2_2l3': lambda: l2 ** 2 * l3,
        's3pil3_2': lambda: s3 * pi * l3 ** 2, 's3pil2_2': lambda: s3 * pi * l2 ** 2, 's3pil2l3': lambda: s3 * pi * l2 * l3,
        'pil2_2': lambda: pi * l2 ** 2, 'pil3_2': lambda: pi * l3 ** 2, 'pil2l3': lambda: pi * l2 * l3,
        'G': lambda: G, 'Gl2': lambda: G * l2, 'Gl3': lambda: G * l3, 'pi3x': lambda: pi ** 3,
        'Lm3_2': lambda: L(chim3, 2), 's3Lm3_2': lambda: s3 * L(chim3, 2),
        'Lm3_2l3': lambda: L(chim3, 2) * l3, 's3Lm3_2l3': lambda: s3 * L(chim3, 2) * l3,
        'L12_3': lambda: L(chi12, 3), 's3L12_3': lambda: s3 * L(chi12, 3),
        'L12_2': lambda: L(chi12, 2), 's3L12_2': lambda: s3 * L(chi12, 2),
        'L5_2': lambda: L(chi5, 2), 'L5_3': lambda: L(chi5, 3), 'Lm8_2': lambda: L(chim8, 2), 'L8_3': lambda: L(chi8, 3),
        'Lm4_4': lambda: L(chim4, 4), 'Lm3_4': lambda: L(chim3, 4), 's3Lm3_4': lambda: s3 * L(chim3, 4),
    }
    return [table[nm]() for nm in names]


def main():
    s, b, n, dps = (int(x) for x in sys.argv[1:5])
    names = sys.argv[5].split(',') if len(sys.argv) > 5 and sys.argv[5] else []
    t0 = time.time()
    T = qvf.bbp_basis(s, b, n, dps)
    tv = target_values(names, dps)
    vals = tv + T
    digits = dps - 20
    rels, dt = qvf.lll_relations(vals, digits)
    goods = qvf.good_relations(rels, digits)
    k = len(names)
    rank, rows = qvf.proj_rank(goods, k) if k else (0, [])
    print(f"s={s} b={b} n={n} dps={dps} targets={names}")
    print(f"relations found: {len(goods)} (lattice dim {len(vals)}), LLL {dt:.1f}s, total {time.time()-t0:.1f}s")
    print(f"rank of target projection: {rank}")
    for r in rows:
        print("  target combo in span:", dict((nm, str(x)) for nm, x in zip(names, r) if x != 0))
    # print relations that involve targets, smallest first
    inv = [g for g in goods if any(g[:k])]
    inv.sort(key=lambda g: max(abs(x) for x in g))
    for g in inv[:6]:
        print("  rel:", g[:k], "|", g[k:])
    # noise check: smallest non-relation residual
    bad = [lr for c, lr, mx in rels if lr >= -(digits - 60)]
    if bad:
        print(f"  smallest non-relation log10 residual: {min(bad):.1f}")


if __name__ == '__main__':
    main()
