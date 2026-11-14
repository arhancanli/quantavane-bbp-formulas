"""Find a sparse certificate: express the target as a rational combination of few five-term relations.
Rebuilds the relation list exactly as sunit5_solve.py does, then tries many column orders in an exact
fraction-based elimination (basic solutions use at most rank-many relations; we keep the sparsest found).
Usage: python -I sparse_cert.py exc.json target_spec out.json [tries]"""
import sys, json, random
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import sunit5_solve as S


def build(exc_path, target_spec):
    exc = json.load(open(exc_path))
    pts = {}
    for e, (re, im) in exc:
        z = complex(re, im)
        pts[S.key(z)] = (tuple(e), z)
    maps = [(lambda z: z, 1), (lambda z: 1 - 1 / z, 1), (lambda z: 1 / (1 - z), 1),
            (lambda z: 1 / z, -1), (lambda z: 1 - z, -1), (lambda z: z / (z - 1), -1)]
    def canon(z):
        if abs(z.imag) < 1e-10:
            return (None, 0)
        best = None
        for f, s in maps:
            v = f(z)
            for w, s2 in ((v, s), (v.conjugate(), -s)):
                k = S.key(w)
                if best is None or k < best[0]:
                    best = (k, s2)
        if best[0] not in pts:
            return None
        return best
    zs = [v[1] for v in pts.values()]
    rels, seen = [], set()
    for x in zs:
        if canon(x) is None:
            continue
        for y in zs:
            if abs(1 - x * y) < 1e-9:
                continue
            w4 = 1 - x * y
            if S.key(w4) not in pts and S.key(1 - w4) not in pts:
                continue
            vec, ok = {}, True
            for a in (x, y, (1 - x) / (1 - x * y), w4, (1 - y) / (1 - x * y)):
                c = canon(a)
                if c is None:
                    ok = False
                    break
                if c[0] is not None:
                    vec[c[0]] = vec.get(c[0], 0) + c[1]
            if not ok:
                continue
            vec = {k: v for k, v in vec.items() if v}
            sig = tuple(sorted(vec.items()))
            if vec and sig not in seen:
                seen.add(sig)
                rels.append((vec, (x, y)))
    target = {}
    for part in target_spec.split(';'):
        pt, c = part.split(':')
        re, im = (float(t) for t in pt.split(','))
        cz = canon(complex(re, im))
        if cz[0] is not None:
            target[cz[0]] = target.get(cz[0], Fr(0)) + Fr(c) * cz[1]
    return rels, target


def solve_order(rels, target, order):
    """Exact elimination on columns = relations (in given order); returns solution dict or None."""
    keys = sorted(set(k for v, _ in rels for k in v) | set(target))
    rows = {k: {} for k in keys}            # row per class: col -> coef
    for j in order:
        for k, v in rels[j][0].items():
            rows[k][j] = Fr(v)
    R = [(dict(rows[k]), target.get(k, Fr(0))) for k in keys]
    piv = []
    used_rows = []
    for j in order:
        r = next((i for i, (row, rhs) in enumerate(R) if i not in used_rows and row.get(j, 0) != 0), None)
        if r is None:
            continue
        row, rhs = R[r]
        pv = row[j]
        row = {c: v / pv for c, v in row.items()}
        rhs = rhs / pv
        R[r] = (row, rhs)
        for i, (row2, rhs2) in enumerate(R):
            if i != r and row2.get(j, 0) != 0:
                f = row2[j]
                nr = dict(row2)
                for c, v in row.items():
                    nv = nr.get(c, 0) - f * v
                    if nv == 0:
                        nr.pop(c, None)
                    else:
                        nr[c] = nv
                R[i] = (nr, rhs2 - f * rhs)
        used_rows.append(r)
        piv.append((j, r))
    for i, (row, rhs) in enumerate(R):
        if i not in used_rows and rhs != 0:
            return None
    sol = {}
    for j, r in piv:
        row, rhs = R[r]
        if rhs != 0:
            sol[j] = rhs                       # free variables set to 0
    return sol


def main():
    rels, target = build(sys.argv[1], sys.argv[2])
    tries = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    print("relations:", len(rels), flush=True)
    best = None
    rng = random.Random(1)
    base = list(range(len(rels)))
    base.sort(key=lambda j: len(rels[j][0]))
    for t in range(tries):
        order = base[:] if t == 0 else rng.sample(base, len(base))
        sol = solve_order(rels, target, order)
        if sol is not None and (best is None or len(sol) < len(best)):
            best = sol
            print(f"try {t}: support {len(sol)}", flush=True)
    out = []
    for j, c in best.items():
        x, y = rels[j][1]
        out.append({"coef": str(c), "x": [x.real, x.imag], "y": [y.real, y.imag]})
    json.dump(out, open(sys.argv[3], 'w'), indent=1)
    print("best support:", len(best))


if __name__ == '__main__':
    main()
