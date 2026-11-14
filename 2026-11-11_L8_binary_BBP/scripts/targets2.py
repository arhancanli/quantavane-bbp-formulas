"""Extended target constants (degree 2 and 3) for BBP span searches."""
import mpmath as mp

def kron(d, n):
    """Kronecker symbol (d/n) for fundamental discriminant d, n >= 1."""
    from sympy import jacobi_symbol
    from sympy.ntheory import factorint
    if n == 0:
        return 0
    res = 1
    for p, e in factorint(n).items():
        if p == 2:
            if d % 2 == 0:
                return 0
            v = 1 if d % 8 in (1, 7) else -1
        else:
            if d % p == 0:
                return 0
            v = jacobi_symbol(d % p, p)
        res *= v ** e
    return res

def Lval(d, s):
    m = abs(d)
    return mp.fsum(kron(d, r) * mp.zeta(s, mp.mpf(r) / m) for r in range(1, m) if kron(d, r)) / mp.mpf(m) ** s

def values(names, dps):
    mp.mp.dps = dps + 30
    pi, l2, l3 = mp.pi, mp.log(2), mp.log(3)
    ls = mp.log(1 + mp.sqrt(2))      # log(1+sqrt2)
    lg = mp.log((1 + mp.sqrt(5)) / 2)  # log golden ratio
    lt = mp.log(2 + mp.sqrt(3))
    out = []
    for nm in names:
        # syntax: [sK*]base ; sK means times sqrt(K)
        mult = mp.mpf(1)
        core = nm
        if nm.startswith('s') and '*' in nm:
            k, core = nm[1:].split('*', 1)
            mult = mp.sqrt(int(k))
        if core.startswith('L(') :
            d, s = core[2:-1].split(',')
            v = Lval(int(d), int(s))
        else:
            env = dict(pi=pi, l2=l2, l3=l3, ls=ls, lg=lg, lt=lt, z3=mp.zeta(3), z5=mp.zeta(5), G=mp.catalan,
                       l5=mp.log(5), l7=mp.log(7))
            v = eval(core, {"__builtins__": {}}, env)
        out.append(mult * v)
    return out
