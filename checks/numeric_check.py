# Independent numerical check of every identity in the paper (written for papers/bbp-formulas-2026; it uses only the
# statements, not the lane's scripts). mpmath; L-values via Hurwitz zeta; polylogarithms via mpmath.polylog.
from mpmath import mp, mpf, mpc, sqrt, pi, log, zeta, polylog, exp, im, arg, atan, ceil, nstr, fabs, clsin
mp.dps = 1100
def L8(): return (zeta(2, mpf(1)/8) + zeta(2, mpf(3)/8) - zeta(2, mpf(5)/8) - zeta(2, mpf(7)/8)) / 64
def L20():
    chi = {1: 1, 3: 1, 7: 1, 9: 1, 11: -1, 13: -1, 17: -1, 19: -1}
    return sum(s * zeta(2, mpf(a) / 20) for a, s in chi.items()) / 400
def Ti2(x): return im(polylog(2, mpc(0, x)))
def D(z): return im(polylog(2, z)) + arg(1 - z) * log(abs(z))
out = []
def rep(name, val): out.append(f'{name}: {nstr(val, 5)}')
l8, l20 = L8(), L20()
out.append('L(2,chi_-8)  = ' + nstr(l8, 40))
out.append('L(2,chi_-20) = ' + nstr(l20, 40))
out.append('Ti_2(1/sqrt2) = ' + nstr(Ti2(1/sqrt(2)), 30))
# Theorem 1.1 (Q(sqrt-2)): 6 Ti2(1/sqrt2) + Ti2(1/(2 sqrt2)) = 4 sqrt2 L8 - (3 pi/4) log 2
rep('|Thm 1 (Ti_2 form)|', 6 * Ti2(1/sqrt(2)) + Ti2(1/(2*sqrt(2))) - (4*sqrt(2)*l8 - 3*pi/4*log(2)))
z8 = exp(mpc(0, 1) * pi / 4)
rep('|Bloch-Wigner form, Q(zeta_8)|', 6*D(mpc(0, 1)/sqrt(2)) + D(mpc(0, 1)/(2*sqrt(2))) - 4*(D(z8) + D(z8**3)))
rep('|D(z8)+D(z8^3) - sqrt2 L8|', D(z8) + D(z8**3) - sqrt(2)*l8)
# BBP form, base 64, length 12
A8 = [32, 0, 8, 0, 8, 0, -4, 0, -1, 0, -1, 0]
P = sum(mpf(64)**(-k) * sum(mpf(A8[j]) / (12*k + j + 1)**2 for j in range(12)) for k in range(0, 600))
rep('|L8 - 3 sqrt2/32 pi log2 - (3/128) P(2,64,12,A)|', l8 - 3*sqrt(2)/32*pi*log(2) - mpf(3)/128*P)
e = [2, 1, 2, -2, -1, -2]
S = sum(mpf(e[n % 6]) / (mpf(2)**n * (2*n + 1)**2) for n in range(0, 3700))
rep('|L8 - 3 sqrt2/32 pi log2 - (3/8) sum e_n/(2^n (2n+1)^2)|', l8 - 3*sqrt(2)/32*pi*log(2) - mpf(3)/8*S)
Sa = sum((-1)**k * (mpf(6)/2**k + mpf(1)/2**(3*k + 1)) / (2*k + 1)**2 for k in range(0, 3700))
rep('|Cor (a)|', Sa - (8*l8 - 3*sqrt(2)/4*pi*log(2)))
# Q(sqrt-5)
z = lambda a: exp(mpc(0, 1) * pi * a / 20) / sqrt(2)
ImL = lambda a: im(polylog(2, z(a)))
rep('|Thm A (Im Li_2 form)|', ImL(1) + ImL(7) + ImL(9) - ImL(17) - (3*sqrt(5)/4*l20 - pi/4*log(2)))
z20 = lambda b: exp(mpc(0, 1) * pi * b / 10)
rep('|Thm A (D form)|', D(z(1)) + D(z(7)) + D(z(9)) - D(z(17)) - mpf(3)/4*sum(D(z20(b)) for b in (1, 3, 7, 9)))
rep('|sum_b D(zeta_20^b) - sqrt5 L20|', sum(D(z20(b)) for b in (1, 3, 7, 9)) - sqrt(5)*l20)
c20 = [1, 1, -1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, -1, 1, 1, 0]
def c(k):
    r = (k - 1) % 40
    return c20[r] if r < 20 else -c20[r - 20]
S20 = sum(mpf(c(k)) / (mpf(2)**((k + 1)//2) * k * k) for k in range(1, 7400))
rep('|Cor A: L20 - sqrt5/15 pi log2 - (4/3) sum c_k/(2^ceil(k/2) k^2)|', l20 - sqrt(5)/15*pi*log(2) - mpf(4)/3*S20)
A = [2**19, 2**19, -2**18, 0, 0, 2**17, 2**16, 0, 2**15, 0, 2**14, 0, 2**13, 2**13, 0, 0, -2**11, 2**11, 2**10, 0,
     -2**9, -2**9, 2**8, 0, 0, -2**7, -2**6, 0, -2**5, 0, -2**4, 0, -2**3, -2**3, 0, 0, 2, -2, -1, 0]
P20 = sum(mpf(2)**(-20*k) * sum(mpf(A[j]) / (40*k + j + 1)**2 for j in range(40)) for k in range(0, 200))
rep('|Cor A: (1/(3 2^18)) P(2,2^20,40,A) form|', l20 - sqrt(5)/15*pi*log(2) - P20/(3*mpf(2)**18))
S1 = sum(mpf(c(k)) / (mpf(2)**((k + 1)//2) * k) for k in range(1, 7400))
rep('|degree 1: sum c_k/(2^ceil(k/2) k) - pi sqrt5/10|', S1 - pi*sqrt(5)/10)
C2 = lambda t: clsin(2, t)
rep('|B1|', ImL(7) + ImL(17) - ((21*C2(2*pi/5) - C2(4*pi/5))/16 - 3*pi/20*log(2)))
rep('|B2|', ImL(1) - ImL(17) - ((6*sqrt(5)*l20 - 11*C2(2*pi/5) - 10*C2(4*pi/5))/16 - pi/40*log(2)))
rep('|B3|', ImL(9) - ImL(17) - ((6*sqrt(5)*l20 - 10*C2(2*pi/5) + 11*C2(4*pi/5))/16 - 3*pi/40*log(2)))
# red controls
rep('red: Thm 1 with 6 -> 5', 5 * Ti2(1/sqrt(2)) + Ti2(1/(2*sqrt(2))) - (4*sqrt(2)*l8 - 3*pi/4*log(2)))
rep('red: Cor A with companion sqrt5/16', l20 - sqrt(5)/16*pi*log(2) - mpf(4)/3*S20)
out.append('vol(H^3/PSL2(Z[sqrt-2])) = (2 sqrt2/3) L8 = ' + nstr(2*sqrt(2)/3*l8, 20))
print('\n'.join(out))
