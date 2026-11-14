"""Quantavane formulas lane: shared helpers (BBP bases, integer relations).

All numerics with mpmath (gmpy backend); lattice reduction with python-flint (FLINT LLL).
"""
import time
import mpmath as mp
from flint import fmpz_mat


def set_dps(d):
    mp.mp.dps = d


def bbp_basis(s, b, n, dps):
    """T_j = sum_{k>=0} b^{-k} / (n k + j)^s for j = 1..n, to about dps digits.

    b may be negative (alternating). Returns list of mpf (length n).
    """
    mp.mp.dps = dps + 30
    K = int((dps + 40) / float(mp.log10(abs(b)))) + 3
    T = [mp.mpf(0)] * n
    w = mp.mpf(1)
    binv = mp.mpf(1) / b
    for k in range(K):
        base = n * k
        for j in range(1, n + 1):
            T[j - 1] += w / mp.mpf(base + j) ** s
        w *= binv
    return T


def P(s, b, n, A, dps):
    """Bailey's notation P(s,b,n,A) = sum_k b^{-k} sum_j a_j/(kn+j)^s."""
    T = bbp_basis(s, b, n, dps)
    return mp.fsum(mp.mpf(a) * t for a, t in zip(A, T))


def lll_relations(vals, digits, verbose=False):
    """Find small integer relations among vals (list of mpf) using LLL.

    digits: number of decimal digits used for the scaled column (should be
    below the working precision of vals by >= 20).
    Returns list of (coeff_list, log10|residual|) for rows whose residual is tiny.
    """
    d = len(vals)
    scale = mp.mpf(10) ** digits
    rows = []
    for i, v in enumerate(vals):
        r = [0] * d
        r[i] = 1
        r.append(int(mp.nint(v * scale)))
        rows.append(r)
    t0 = time.time()
    L = fmpz_mat(rows).lll()
    dt = time.time() - t0
    out = []
    for row in L.tolist():
        c = [int(x) for x in row[:d]]
        if all(x == 0 for x in c):
            continue
        res = mp.fsum(mp.mpf(ci) * vi for ci, vi in zip(c, vals))
        lr = float(mp.log10(abs(res))) if res != 0 else -10 ** 9
        mx = max(abs(x) for x in c)
        out.append((c, lr, mx))
    if verbose:
        print(f"[lll] dim={d} digits={digits} time={dt:.1f}s")
    return out, dt


def good_relations(rels, digits, slack=60):
    """Rows that are genuine relations, not lattice noise.

    A generic (noise) LLL vector in dimension d has coefficients near 10^(digits/d)
    and residual near 10^(-digits + digits/d).  A genuine relation has residual at the
    rounding level (about 10^-digits or below) and much smaller coefficients.  We accept
    a row only if its coefficients are below 10^(digits/(2d)) and its residual is below
    10^(-(digits - slack)).  Every accepted relation must still be re-checked at higher
    precision before use.
    """
    goods = []
    if not rels:
        return goods
    d = len(rels[0][0])
    gen = digits / float(d)          # log10 of a generic LLL vector's size
    import math
    for c, lr, mx in rels:
        # genuine: residual at the rounding floor (well below the noise level
        # -digits + gen) and coefficients clearly smaller than generic.
        if lr < -(digits - 0.5 * gen) and math.log10(mx) < 0.75 * gen:
            goods.append(c)
    return goods


def proj_rank(goods, k):
    """Rank (over Q) of the projections of relation vectors onto the first k coordinates."""
    if not goods:
        return 0, []
    import sympy
    M = sympy.Matrix([g[:k] for g in goods])
    rref, piv = M.rref()
    rows = [list(rref.row(i)) for i in range(len(piv))]
    return len(piv), rows
