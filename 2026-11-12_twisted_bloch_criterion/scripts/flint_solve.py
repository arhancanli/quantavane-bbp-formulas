"""Exact solve of the OMP-selected system with FLINT: find a sparse rational combination of the selected five-term
relations whose reduced vector equals the target; write it in the format of verify_cert_exact.py.
Usage: python -I flint_solve.py system.omp.json out_cert.json [tries]"""
import sys, json, random
from fractions import Fraction as Fr
from flint import fmpq_mat, fmpq

D = json.load(open(sys.argv[1]))
tries = int(sys.argv[3]) if len(sys.argv) > 3 else 200
S = D['S']
cols = [{(k0, k1): v for k0, k1, v in rv} for rv in D['RV']]
target = {(k0, k1): Fr(v) for k0, k1, v in D['target']}
keys = sorted(set(k for c in cols for k in c) | set(target))
kid = {k: i for i, k in enumerate(keys)}
m, n = len(keys), len(cols)
print(f"system: {m} classes x {n} relations", flush=True)


def solve(order):
    """rref of [A_order | b]; solution with free variables 0. Returns dict col -> Fraction or None."""
    A = fmpq_mat(m, n + 1)
    for jj, j in enumerate(order):
        for k, v in cols[j].items():
            A[kid[k], jj] = v
    for k, v in target.items():
        A[kid[k], n] = fmpq(v.numerator, v.denominator)
    R, rk = A.rref()
    sol = {}
    for i in range(rk):
        piv = next(c for c in range(n + 1) if R[i, c] != 0)
        if piv == n:
            return None            # inconsistent
        val = R[i, n]
        if val != 0:
            sol[order[piv]] = Fr(int(val.p), int(val.q))
    return sol


rng = random.Random(11)
best = None
base = sorted(range(n), key=lambda j: (len(cols[j]), j))
for t in range(tries):
    order = base if t == 0 else rng.sample(range(n), n)
    sol = solve(order)
    if sol is None:
        print("inconsistent: target not in span"); sys.exit(1)
    if best is None or len(sol) < len(best):
        best = sol
        print(f"try {t}: support {len(sol)}", flush=True)
# greedy pruning with the fast solver: restrict to the support and drop columns while solvable
supp = sorted(best)
improved = True
while improved:
    improved = False
    for j in list(supp):
        trial = [c for c in supp if c != j]
        # solve restricted system
        sub_cols = trial
        A = fmpq_mat(m, len(sub_cols) + 1)
        for jj, c in enumerate(sub_cols):
            for k, v in cols[c].items():
                A[kid[k], jj] = v
        for k, v in target.items():
            A[kid[k], len(sub_cols)] = fmpq(v.numerator, v.denominator)
        R, rk = A.rref()
        ok = True; sol = {}
        for i in range(rk):
            piv = next(c for c in range(len(sub_cols) + 1) if R[i, c] != 0)
            if piv == len(sub_cols):
                ok = False; break
            if R[i, len(sub_cols)] != 0:
                v = R[i, len(sub_cols)]
                sol[sub_cols[piv]] = Fr(int(v.p), int(v.q))
        if ok:
            supp = sorted(sol); best = sol; improved = True
            print(f"pruned to {len(supp)}", flush=True)
            break
# exact class-level check
tot = {}
for j, c in best.items():
    for k, v in cols[j].items():
        tot[k] = tot.get(k, 0) + c * v
tot = {k: v for k, v in tot.items() if v}
print("class-level check:", "OK" if tot == target else "MISMATCH", flush=True)
json.dump({'M': D['M'], 'cert': [{'coef': str(c), 'x': D['rels'][j][0], 'y': D['rels'][j][1]} for j, c in best.items()],
           'recipes': D['recipes'], 'target': D['target_cfg']}, open(sys.argv[2], 'w'))
print(f"written {sys.argv[2]} with {len(best)} relations")
