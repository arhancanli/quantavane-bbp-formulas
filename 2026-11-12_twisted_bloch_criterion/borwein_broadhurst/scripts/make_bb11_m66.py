"""Configs for Borwein-Broadhurst (28) D=-11, (29) D=-15, (37) D=-24 (and (30) D=-20) as five-term certificate targets."""
import json, cmath, math
from fractions import Fraction as Fr
from sympy import jacobi_symbol

def kron(d, n):
    from sympy.ntheory import factorint
    res = 1
    for p, e in factorint(n).items():
        if p == 2:
            v = 0 if d % 2 == 0 else (1 if d % 8 in (1, 7) else -1)
        else:
            v = 0 if d % p == 0 else jacobi_symbol(d % p, p)
        res *= v ** e
    return res

def gauss(D, M):
    """sqrt(D) (D<0 fundamental, |D| | M) as a dict exponent -> coefficient in Q(zeta_M): sum_t chi_D(t) zeta_|D|^t."""
    n = -D; out = {}
    for t in range(1, n):
        c = kron(D, t)
        if c:
            k = (M // n) * t % M
            out[k] = out.get(k, 0) + c
    return {k: Fr(v) for k, v in out.items() if v}

def mul(a, b, M):
    r = {}
    for i, x in a.items():
        for j, y in b.items():
            k = (i + j) % M
            r[k] = r.get(k, 0) + x * y
    return {k: v for k, v in r.items() if v}

def num(d, M):
    return sum(complex(float(c)) * cmath.exp(2j * math.pi * k / M) for k, c in d.items())

cases = {
 # D: (M, S, x, y, sqrtD_scale, ks_needed, coeffs dict k->c, Z: list of (exponent of zeta_M, coefficient))
 -11: dict(M=66, S=[2, 3, 11], x=Fr(-5, 6), y=Fr(1, 6), coeffs={1: 15, 2: -10, 5: 1}),
}
for D, c in cases.items():
    M = c['M']
    g = gauss(D, M)                      # sqrt(D) = i sqrt|D|
    sD = num(g, M); assert abs(sD - 1j * math.sqrt(-D)) < 1e-9, (D, sD)
    w = {0: c['x']}
    for k, v in g.items():
        w[k] = w.get(k, 0) + c['y'] * v
    wv = num(w, M); assert abs(abs(wv) - 1) < 1e-12, (D, wv)
    pw = {1: w}
    for k in range(2, max(c['coeffs']) + 1):
        pw[k] = mul(pw[k - 1], w, M)
    # numeric check of the identity: sum coeff_k Cl2(k theta) = |D| sum_{k=1}^{|D|/2} chi_D(k) Cl2(2 pi k/|D|)
    import mpmath as mp
    mp.mp.dps = 30
    th = mp.arg(mp.mpc(wv.real, wv.imag))
    lhs = sum(cf * mp.clsin(2, k * th) for k, cf in c['coeffs'].items())
    n = -D
    rhs = n * sum(kron(D, k) * mp.clsin(2, 2 * mp.pi * k / n) for k in range(1, n // 2 + 1))
    print(D, 'theta', th, 'lhs', lhs, 'rhs', rhs, 'diff', lhs - rhs)
    target = [["pow", {str(k): str(v) for k, v in pw[k].items()}, str(cf)] for k, cf in c['coeffs'].items()]
    for k in range(1, n // 2 + 1):
        ch = kron(D, k)
        if ch:
            target.append([(M // n) * k % M, 0, str(-n * ch)])
    cfg = {"M": M, "S": c['S'], "rounds": 2, "cap": 1200, "ks": [0], "tries": 6,
           "extra": [{str(k): str(v) for k, v in pw[j].items()} for j in sorted(pw)], "target": target}
    json.dump(cfg, open(f'search_config_bb{-D}_m66.json', 'w'))
print('configs written')
