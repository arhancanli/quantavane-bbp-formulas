"""Print a certificate (from cert_search.py) in readable form, and check each five-term relation numerically (50 digits).
Usage: python -I describe_cert.py cert.json"""
import sys, json
import mpmath as mp
mp.mp.dps = 50
d = json.load(open(sys.argv[1])); M = d['M']; R = d['recipes']
SYMN = ['{}', '1 - 1/({})', '1/(1 - {})', '1/({})', '1 - ({})', '({0})/({0} - 1)']
def desc(i):
    r = R[i]
    if r[0] == 'bbp':
        a, k = r[1], r[2]
        s = f"z40^{a}" if a else "1"
        return s + (f"*2^(-{k}/2)" if k else "")
    w = desc(r[1]); s = SYMN[r[2]].format(w)
    return f"conj({s})" if r[3] else s
def val(i):
    r = R[i]
    if r[0] == 'bbp':
        return mp.expjpi(mp.mpf(2 * r[1]) / M) / mp.sqrt(2) ** r[2]
    w = val(r[1]); si = r[2]
    x = [w, 1 - 1 / w, 1 / (1 - w), 1 / w, 1 - w, w / (w - 1)][si]
    return mp.conj(x) if r[3] else x
def D(z): return mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z))
worst = 0
for it in d['cert']:
    x, y = val(it['x']), val(it['y'])
    q = 1 - x * y
    res = D(x) + D(y) + D((1 - x) / q) + D(q) + D((1 - y) / q)
    worst = max(worst, abs(res))
    print(f"  {it['coef']:>8}  x = {desc(it['x'])}   y = {desc(it['y'])}")
print(f"max |five-term residual| over the certificate (50 digits): {mp.nstr(worst, 3)}")
