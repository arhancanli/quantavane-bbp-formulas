"""Is a target D-combination in the Q-span of found five-term relations plus distribution relations?
Usage: python -I solve_d.py rels.json "r,kk:coef;r,kk:coef;..."   (N = 24, keys as in fiveterm.py)"""
import sys, json, cmath, math
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fractions import Fraction as Fr
import sympy
from fiveterm import canon, N, KMAX


def dist_relations():
    """D(z^n) = n * sum_{w^n = 1} D(w z), for z = zeta^r 2^(-kk/2), n = 2, 3."""
    rels = []
    for r in range(N):
        for kk in range(0, KMAX + 1):
            z = cmath.exp(2j * math.pi * r / N) * 2 ** (-kk / 2)
            for n in (2, 3):
                if n * kk > KMAX:
                    continue
                vec = {}
                c = canon(z ** n)
                if c is None:
                    continue
                if c[0] is not None:
                    vec[c[0]] = vec.get(c[0], 0) + c[1]
                ok = True
                for j in range(n):
                    w = cmath.exp(2j * math.pi * j / n) * z
                    cw = canon(w)
                    if cw is None:
                        ok = False
                        break
                    if cw[0] is not None:
                        vec[cw[0]] = vec.get(cw[0], 0) - n * cw[1]
                if ok:
                    vec = {k: v for k, v in vec.items() if v}
                    if vec:
                        rels.append(vec)
    return rels


def main():
    data = json.load(open(sys.argv[1]))
    rels = [{tuple(k): v for k, v in d["rel"]} for d in data]
    nfive = len(rels)
    rels += dist_relations()
    target = {}
    for part in sys.argv[2].split(';'):
        key, c = part.split(':')
        r, kk = (int(t) for t in key.split(','))
        target[(r, kk)] = Fr(c)
    keys = sorted(set(k for rel in rels for k in rel) | set(target))
    idx = {k: i for i, k in enumerate(keys)}
    M = sympy.Matrix([[rel.get(k, 0) for k in keys] for rel in rels])
    t = sympy.Matrix([[sympy.Rational(target.get(k, 0).numerator, target.get(k, 0).denominator) for k in keys]])
    print(f"relations: {nfive} five-term + {len(rels)-nfive} distribution; unknowns {len(keys)}; rank {M.rank()}")
    # solve c^T M = t
    A = M.T
    try:
        sol, params = A.gauss_jordan_solve(t.T)
        print("TARGET IS IN THE SPAN")
        sol = sol.subs({p: 0 for p in params})
        used = [(i, sol[i]) for i in range(len(rels)) if sol[i] != 0]
        for i, c in used:
            kind = 'five-term' if i < nfive else 'distribution'
            extra = f" x={data[i]['x']} y={data[i]['y']}" if i < nfive else ''
            print(f"  {c} * [{kind} #{i}] {rels[i]}{extra}")
    except ValueError:
        print("target NOT in span")


if __name__ == '__main__':
    main()
