"""Five-term-relation certificate search in Q(zeta_M) (numeric search, exact recipes for later verification).

Elements are stored by their vector of embeddings sigma_t(w), t in (Z/M)^*, plus a recipe from which the exact element
can be rebuilt.  Exceptional S-units (w and 1-w with norms supported on S) are closed under the 12 Bloch-Wigner
symmetries; five-term relations R(x,y) are collected when 1-xy is an S-unit; finally an exact rational combination of
relations equal to the target D-combination is sought (exact sparse elimination with Fractions).
Usage: python -I cert_search.py config.json out.json
"""
import sys, json, math, cmath, time, random
from math import gcd
from fractions import Fraction as Fr
import numpy as np

cfg = json.load(open(sys.argv[1]))
M = cfg['M']; S = cfg['S']; BOX = cfg.get('box', 16); ROUNDS = cfg.get('rounds', 1); CAP = cfg.get('cap', 4000)
U = [t for t in range(1, M) if gcd(t, M) == 1]
NU = len(U)
ZE = np.array([[cmath.exp(2j * math.pi * j * t / M) for j in range(M)] for t in U])   # sigma_t(zeta^j)
LOGS = np.array([math.log(p) for p in S])

def emb_from_powers(cs):
    """cs: dict exponent -> coefficient (float ok).  returns embedding vector."""
    v = np.zeros(NU, dtype=complex)
    for j, c in cs.items():
        v += c * ZE[:, j % M]
    return v

def is_Sunit_logn(ln):
    """ln: log|N(w)| (float).  True if it equals sum e_p log p with integer e_p in a box (exact lattice check)."""
    # brute force small box (|e_p| <= 60)
    if len(S) == 1:
        e = ln / LOGS[0]
        return abs(e - round(e)) < 1e-7
    rng = range(-60, 61)
    if len(S) == 2:
        for e2 in rng:
            r = ln - e2 * LOGS[1]
            e1 = r / LOGS[0]
            if abs(e1 - round(e1)) < 1e-7 and abs(e1) <= 400:
                return True
        return False
    if len(S) == 3:
        for e2 in range(-BOX, BOX + 1):
            for e3 in range(-BOX, BOX + 1):
                r = ln - e2 * LOGS[1] - e3 * LOGS[2]
                e1 = r / LOGS[0]
                if abs(e1 - round(e1)) < 1e-7 and abs(e1) <= 400:
                    return True
        return False
    raise ValueError

def logn(v):
    return float(np.sum(np.log(np.abs(v))))

# ---------- element store ----------
vecs = []        # embedding vectors
recipes = []     # recipes
index = {}       # key -> id

def key(v):
    z = v[0]
    return (round(z.real, 8) + 0.0, round(z.imag, 8) + 0.0)

def add(v, rec):
    k = key(v)
    if k in index:
        return index[k], False
    index[k] = len(vecs)
    vecs.append(v); recipes.append(rec)
    return index[k], True

ONE = np.ones(NU, dtype=complex)
SYM = [  # (function on vectors, sign, name)
    (lambda w: w, 1, 'id'), (lambda w: 1 - 1 / w, 1, 'a'), (lambda w: 1 / (1 - w), 1, 'b'),
    (lambda w: 1 / w, -1, 'inv'), (lambda w: 1 - w, -1, 'om'), (lambda w: w / (w - 1), -1, 'c')]
# complex conjugation: sigma_{-t} ; in our U ordering, conj vector = np.conj(v) (CM field)

def canon(v):
    """canonical class key and sign for D(w) at the identity embedding; None for real w."""
    if abs(v[0].imag) < 1e-10:
        return None, 0
    best = None
    for f, s, _ in SYM:
        w = f(v)
        for ww, ss in ((w, s), (np.conj(w), -s)):
            k = key(ww)
            if best is None or k < best[0]:
                best = (k, ss)
    return best

def exceptional(v):
    if np.any(np.abs(v) < 1e-12) or np.any(np.abs(1 - v) < 1e-12):
        return False
    return is_Sunit_logn(logn(v)) and is_Sunit_logn(logn(1 - v))

# ---------- seeds ----------
def seed_list():
    out = []
    s2 = {M // 8: 0.5, (M - M // 8): 0.5}  # 1/sqrt2 = (zeta^(M/8) + zeta^(-M/8))/2  (requires 8 | M)
    for a in range(M):
        for k in cfg.get('ks', [0, 1, 2, 3, 4]):
            # z = zeta^a * (1/sqrt2)^k
            cs = {a: 1.0}
            for _ in range(k):
                new = {}
                for j, c in cs.items():
                    for jj, cc in s2.items():
                        new[(j + jj) % M] = new.get((j + jj) % M, 0) + c * cc
                cs = new
            if k == 0 and a == 0:
                continue
            out.append((emb_from_powers(cs), ('bbp', a, k)))
    for extra in cfg.get('extra', []):
        cs = {int(j): float(Fr(c)) for j, c in extra.items()}
        out.append((emb_from_powers(cs), ('pow', extra)))
    return out

t0 = time.time()
for v, rec in seed_list():
    if exceptional(v):
        i, new = add(v, rec)
print(f"seeds exceptional: {len(vecs)}  ({time.time()-t0:.1f}s)", flush=True)

def close_sym(ids):
    added = []
    for i in ids:
        for si, (f, s, nm) in enumerate(SYM):
            w = f(vecs[i])
            for cj in (0, 1):
                ww = np.conj(w) if cj else w
                j, new = add(ww, ('sym', i, si, cj))
                if new:
                    added.append(j)
    return added

base = list(range(len(vecs)))
close_sym(base)
print(f"after symmetry closure: {len(vecs)}", flush=True)

rels = []      # (x_id, y_id)
relset = set()
def sunit_mask(ln):
    """vectorized: ln array of log|N|; True where ln is in the lattice spanned by log p, p in S (box)."""
    ok = np.zeros(ln.shape, dtype=bool)
    if len(S) == 1:
        e = ln / LOGS[0]
        return np.abs(e - np.round(e)) < 1e-7
    import itertools as it
    rest = LOGS[1:]
    boxes = [range(-BOX, BOX + 1)] * (len(S) - 1)
    for es in it.product(*boxes):
        r = (ln - sum(e * l for e, l in zip(es, rest))) / LOGS[0]
        ok |= np.abs(r - np.round(r)) < 1e-7
    return ok

def find_relations(xs, ys):
    Y = np.array([vecs[j] for j in ys])           # (ny, NU)
    for xi in xs:
        X = vecs[xi]
        P = 1 - Y * X                              # 1 - x y
        ok = np.all(np.abs(P) > 1e-10, axis=1)
        ln = np.sum(np.log(np.abs(np.where(ok[:, None], P, 1))), axis=1)
        good = ok & sunit_mask(ln)
        for r in np.nonzero(good)[0]:
            yi = ys[r]
            kk = (min(xi, yi), max(xi, yi))
            if kk in relset:
                continue
            relset.add(kk)
            rels.append((xi, yi))

for rd in range(ROUNDS):
    ids = list(range(len(vecs)))
    if len(ids) > CAP:
        random.Random(rd).shuffle(ids); ids = sorted(ids[:CAP])
    t1 = time.time()
    n0 = len(rels)
    find_relations(ids, ids)
    print(f"round {rd}: elements {len(vecs)}, relations {len(rels)} (+{len(rels)-n0}) ({time.time()-t1:.1f}s)", flush=True)
    if rd < ROUNDS - 1:
        newids = []
        for (xi, yi) in rels[n0:]:
            x, y = vecs[xi], vecs[yi]
            d = 1 - x * y
            for j, w in ((2, (1 - x) / d), (3, d), (4, (1 - y) / d)):
                k2, new = add(w, ('rel', xi, yi, j))
                if new:
                    newids.append(k2)
        close_sym(newids)
        print(f"   new elements {len(vecs)}", flush=True)

# ---------- linear algebra ----------
def rel_vector(xi, yi):
    x, y = vecs[xi], vecs[yi]
    d = 1 - x * y
    vec = {}
    for w in (x, y, (1 - x) / d, d, (1 - y) / d):
        k, s = canon(w)
        if k is None:
            continue
        vec[k] = vec.get(k, 0) + s
    return {k: v for k, v in vec.items() if v}

target = {}
for item in cfg['target']:
    if item[0] == 'pow':
        cs = {int(j): float(Fr(c)) for j, c in item[1].items()}
        kk, sg = canon(emb_from_powers(cs))
        if kk is not None:
            target[kk] = target.get(kk, Fr(0)) + Fr(item[2]) * sg
        continue
    a, k, c = item  # point zeta^a * (1/sqrt2)^k with coefficient c
    cs = {a: 1.0}
    s2 = {M // 8: 0.5, (M - M // 8): 0.5}
    for _ in range(k):
        new = {}
        for j, cc in cs.items():
            for jj, c2 in s2.items():
                new[(j + jj) % M] = new.get((j + jj) % M, 0) + cc * c2
        cs = new
    kk, s = canon(emb_from_powers(cs))
    target[kk] = target.get(kk, Fr(0)) + Fr(c) * s
target = {k: v for k, v in target.items() if v}
RV = [rel_vector(*r) for r in rels]
keys = sorted(set(k for v in RV for k in v) | set(target))
missing = [k for k in target if k not in set(k2 for v in RV for k2 in v)]
print(f"classes {len(keys)}, relations {len(RV)}; target classes not touched by any relation: {len(missing)}", flush=True)

# numeric feasibility check (least squares)
import scipy.sparse as sp, scipy.sparse.linalg as spl
kid = {k: i for i, k in enumerate(keys)}
rows, cols, data = [], [], []
for j, v in enumerate(RV):
    for k, c in v.items():
        rows.append(kid[k]); cols.append(j); data.append(float(c))
A = sp.csr_matrix((data, (rows, cols)), shape=(len(keys), len(RV)))
b = np.zeros(len(keys))
for k, c in target.items():
    b[kid[k]] = float(c)
sol = spl.lsqr(A, b, atol=1e-14, btol=1e-14, iter_lim=20000)
res = np.linalg.norm(A @ sol[0] - b)
print(f"least-squares residual {res:.3e}", flush=True)
pass

# ---------- exact sparse solve (Fractions) ----------
def exact_solve(order):
    basis = {}   # pivot class -> (vec dict, expr dict)
    for j in order:
        v = {k: Fr(c) for k, c in RV[j].items()}
        e = {j: Fr(1)}
        # reduce
        changed = True
        while v:
            piv = next((k for k in v if k in basis), None)
            if piv is None:
                break
            bv, be = basis[piv]
            f = v[piv] / bv[piv]
            for k, c in bv.items():
                nv = v.get(k, 0) - f * c
                if nv == 0: v.pop(k, None)
                else: v[k] = nv
            for r, c in be.items():
                ne = e.get(r, 0) - f * c
                if ne == 0: e.pop(r, None)
                else: e[r] = ne
        if v:
            pk = min(v)  # deterministic pivot
            basis[pk] = (v, e)
    # express target
    t = dict(target); ex = {}
    while t:
        piv = next((k for k in t if k in basis), None)
        if piv is None:
            return None
        bv, be = basis[piv]
        f = t[piv] / bv[piv]
        for k, c in bv.items():
            nv = t.get(k, 0) - f * c
            if nv == 0: t.pop(k, None)
            else: t[k] = nv
        for r, c in be.items():
            ne = ex.get(r, 0) + f * c
            if ne == 0: ex.pop(r, None)
            else: ex[r] = ne
    return ex

if res < 1e-8:
    # orthogonal matching pursuit to find a small set of relations whose span contains the target
    t3 = time.time()
    Acsc = A.tocsc()
    norms = np.sqrt(np.asarray(Acsc.multiply(Acsc).sum(axis=0)).ravel()); norms[norms == 0] = 1
    r = b.copy(); S = []
    for it in range(1300):
        c = (Acsc.T @ r) / norms
        if S:
            c[S] = 0
        j = int(np.argmax(np.abs(c)))
        S.append(j)
        As = Acsc[:, S].toarray()
        x, *_ = np.linalg.lstsq(As, b, rcond=None)
        r = b - As @ x
        if np.linalg.norm(r) < 1e-13:
            break
    print(f"OMP: support {len(S)}, residual {np.linalg.norm(r):.2e} ({time.time()-t3:.1f}s)", flush=True)
    keys_S = sorted(set(k for j in S for k in RV[j]) | set(target))
    json.dump({'M': M, 'S': S, 'rels': [list(rels[j]) for j in S],
               'RV': [[[str(k[0]), str(k[1]), int(v)] for k, v in RV[j].items()] for j in S],
               'target': [[str(k[0]), str(k[1]), str(v)] for k, v in target.items()],
               'recipes': recipes, 'target_cfg': cfg['target']}, open(sys.argv[2] + '.omp.json', 'w'))
    print("dumped OMP system", flush=True)
    sys.exit(0)
    best = None
    rng = random.Random(7)
    for tr in range(cfg.get('tries', 30)):
        order = sorted(S, key=lambda j: (len(RV[j]), j)) if tr == 0 else rng.sample(S, len(S))
        ex = exact_solve(order)
        if ex is not None and (best is None or len(ex) < len(best)):
            best = ex
            print(f"  try {tr}: certificate support {len(ex)}", flush=True)
    if best is None:
        print("exact solve on the OMP support failed"); sys.exit(1)
    # prune: drop relations one at a time while the target stays representable
    changed = True
    while changed:
        changed = False
        for j in list(best):
            trial = [k for k in best if k != j]
            ex = exact_solve(trial)
            if ex is not None:
                best = ex; changed = True
                print(f"  pruned to {len(best)}", flush=True)
                break
    tot = {}
    for j, c in best.items():
        for k, v in RV[j].items():
            tot[k] = tot.get(k, 0) + c * v
    tot = {k: v for k, v in tot.items() if v}
    print("class-level check:", "OK" if tot == target else "MISMATCH", flush=True)
    json.dump({'M': M, 'cert': [{'coef': str(c), 'x': rels[j][0], 'y': rels[j][1]} for j, c in best.items()],
               'recipes': recipes, 'target': cfg['target']}, open(sys.argv[2], 'w'))
    print("written", sys.argv[2])
