"""D-value relation search in Q(zeta_M): find all integer relations among Bloch-Wigner values D(z_i) at BBP points
z = zeta_M^a * b^(-k/2) and a basis of Clausen values Cl_2(2 pi a / M).
Usage: python -I dsearch.py M b "a1:k1,a2:k2,..." digits"""
import sys, time
import mpmath as mp
from flint import fmpz_mat

def D(z):
    return mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z))

def lll(vals, digits):
    d = len(vals)
    sc = mp.mpf(10) ** digits
    rows = []
    for i, v in enumerate(vals):
        r = [0] * d; r[i] = 1; r.append(int(mp.nint(v * sc))); rows.append(r)
    L = fmpz_mat(rows).lll()
    out = []
    for row in L.tolist():
        c = [int(x) for x in row[:d]]
        res = mp.fsum(ci * vi for ci, vi in zip(c, vals))
        lr = float(mp.log10(abs(res))) if res != 0 else -1e9
        mx = max(abs(x) for x in c)
        if mx > 0 and lr < -(digits - 15) and mx < 10 ** (digits / (2.0 * d)):
            out.append(c)
    return out

def clausen_basis(M, digits):
    vals = [mp.clsin(2, 2 * mp.pi * a / M) for a in range(1, (M + 1) // 2)]
    idx = list(range(1, (M + 1) // 2))
    rels = lll(vals, digits)
    # greedy: drop variables until no relations remain (exact rank via sympy)
    import sympy
    R = sympy.Matrix(rels) if rels else sympy.zeros(0, len(vals))
    rank = R.rank() if rels else 0
    # choose basis = non-pivot columns of rref of relations (pivot columns are dependent)
    if rels:
        rref, piv = R.rref()
    else:
        piv = ()
    basis = [i for i in range(len(vals)) if i not in piv]
    return [idx[i] for i in basis], [vals[i] for i in basis], len(vals) - rank

def main():
    M = int(sys.argv[1]); b = int(sys.argv[2])
    pts = [tuple(int(x) for x in p.split(':')) for p in sys.argv[3].split(',')]
    digits = int(sys.argv[4])
    mp.mp.dps = digits + 30
    t0 = time.time()
    ca, cv, dim = clausen_basis(M, digits)
    print(f"M={M}: Clausen basis a={ca} (dim {dim})", flush=True)
    sb = mp.sqrt(b)
    dv = [D(mp.expjpi(mp.mpf(2 * a) / M) / sb ** k) for a, k in pts]
    vals = dv + cv
    rels = lll(vals, digits)
    import sympy
    np_ = len(pts)
    print(f"relations: {len(rels)} ({time.time()-t0:.1f}s)")
    if rels:
        R = sympy.Matrix(rels)
        rref, piv = R.rref()
        for i in range(len(piv)):
            row = list(rref.row(i))
            s = " + ".join(f"({row[j]})D[{pts[j][0]}:{pts[j][1]}]" for j in range(np_) if row[j] != 0)
            s2 = " + ".join(f"({row[np_+j]})Cl[{ca[j]}]" for j in range(len(ca)) if row[np_ + j] != 0)
            print("  ", s, "|", s2)

if __name__ == '__main__':
    main()
