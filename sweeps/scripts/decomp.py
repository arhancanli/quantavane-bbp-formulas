"""Decompose a BBP sum P(s, b, n, A) into polylogarithms at the points zeta_n^r * b^(-1/n):
P = sum_r gamma_r Li_s(zeta_n^r x),  x = b^(-1/n),  gamma_r = (1/n) sum_j a_j b^(j/n) zeta_n^(-r j).
Prints gamma_r numerically (and checks the reconstruction).  Usage: python -I decomp.py s b n "a1,...,an" [c0]"""
import sys
import mpmath as mp

def main():
    s, b, n = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    A = [int(t) for t in sys.argv[4].split(',')]
    c0 = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    mp.mp.dps = 60
    x = mp.mpf(b) ** (mp.mpf(-1) / n)
    zeta = mp.expjpi(mp.mpf(2) / n)
    gam = []
    for r in range(n):
        g = mp.fsum(A[j - 1] * mp.mpf(b) ** (mp.mpf(j) / n) * zeta ** (-r * j) for j in range(1, n + 1)) / n / c0
        gam.append(g)
    val = mp.fsum(g * mp.polylog(s, zeta ** r * x) for r, g in enumerate(gam))
    direct = mp.fsum(mp.mpf(A[j-1]) * mp.nsum(lambda k: mp.mpf(b) ** (-k) / (n * k + j) ** s, [0, mp.inf]) for j in range(1, n + 1)) / c0
    print("reconstruction error:", mp.nstr(abs(val - direct), 5), " value:", mp.nstr(direct.real if hasattr(direct,'real') else direct, 30))
    for r, g in enumerate(gam):
        if abs(g) > mp.mpf(10) ** -40:
            print(f"r={r:3d}  angle={r}*2pi/{n}  gamma = {mp.nstr(g.real, 20)} + {mp.nstr(g.imag, 20)} i")

if __name__ == '__main__':
    main()
