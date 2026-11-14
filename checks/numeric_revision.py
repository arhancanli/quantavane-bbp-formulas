# Numerical checks (not proofs) of the statements added in the revision of 12 Nov 2026: Corollary 3.2 (Clausen form),
# the null identity and its equivalence with Almkvist's identity 1.2, the formula for Ti_2(tan t), the real ladder
# behind Broadhurst's eq. (16), Proposition 5.1 (vectors V_1, V_2), the factor 1/2 in Borwein-Broadhurst (27) versus
# BBBZ (10), the degree-one formulas with the same vectors, and the two combinations T_1, T_2 used in Section 5.
# Written from the statements in the paper; mpmath only.  Run: python -I numeric_revision.py
from mpmath import mp, mpf, mpc, sqrt, pi, log, zeta, polylog, exp, im, atan, clsin, nstr
mp.dps = 620
Cl2 = lambda t: clsin(2, t)
L8 = (zeta(2, mpf(1)/8) + zeta(2, mpf(3)/8) - zeta(2, mpf(5)/8) - zeta(2, mpf(7)/8)) / 64
th2 = 2 * atan(sqrt(2))
phi = atan(1 / sqrt(2))
Ti2 = lambda x: im(polylog(2, mpc(0, x)))
out = []
def rep(name, val, digits=5):
    out.append(f'{name}: {nstr(val, digits)}')
mp.dps = 220
rep('Cor 3.2: 6Cl2(t)-(5/2)Cl2(2t)+(1/4)Cl2(4t) - 4 sqrt2 L8', 6*Cl2(th2) - mpf(5)/2*Cl2(2*th2) + Cl2(4*th2)/4 - 4*sqrt(2)*L8)
rep('null: 3Cl2(t)+Cl2(2t)+Cl2(3t)-Cl2(4t)', 3*Cl2(th2) + Cl2(2*th2) + Cl2(3*th2) - Cl2(4*th2))
rep('Almkvist 1.2 with tan(a)=1/sqrt2', Cl2(6*phi - pi) + Cl2(pi + 2*phi) - 2*Cl2(2*phi) + 2*Cl2(pi - 4*phi))
rep('pi - 2 phi - theta_2', pi - 2*phi - th2)
for t in (phi, pi/2 - 2*phi):
    rep(f'Ti2(tan t) - t log tan t - (Cl2(2t)+Cl2(pi-2t))/2 at t = {nstr(t, 8)}',
        Ti2(mp.tan(t)) - t*log(mp.tan(t)) - (Cl2(2*t) + Cl2(pi - 2*t))/2)
rep('duplication Cl2(2x) - 2Cl2(x) + 2Cl2(pi-x) at x = theta_2', Cl2(2*th2) - 2*Cl2(th2) + 2*Cl2(pi - th2))
rep('real ladder 6Li2(-1/2) - Li2(-1/8) + pi^2/3 - (3/2) log^2 2', 6*polylog(2, mpf(-1)/2) - polylog(2, mpf(-1)/8) + pi**2/3 - mpf(3)/2*log(2)**2)
h = mpc(0, 1) / sqrt(2)
rep('Broadhurst (16): Re(6Li2(h)-Li2(h^3)) + pi^2/12 - (3/8)log^2 2', (6*polylog(2, h) - polylog(2, h**3)).real + pi**2/12 - mpf(3)/8*log(2)**2)
lhs27 = 27*Cl2(th2) - 9*Cl2(2*th2) + Cl2(3*th2)
rhs27 = 8*(Cl2(pi/4) + Cl2(3*pi/4))
rep('BB98 (27): (1/2)(27Cl2-9Cl2+Cl2) - 8(Cl2(pi/4)+Cl2(3pi/4))', lhs27/2 - rhs27)
rep('BBBZ (10) as printed: LHS/RHS (should be 2 if the 1/2 is missing)', lhs27/rhs27, 20)
# degree-one formulas with the same vectors
def P(s, b, n, A, terms):
    return sum(mpf(b)**(-q) * sum(mpf(A[j]) / (n*q + j + 1)**s for j in range(n)) for q in range(terms))
c20 = [1, 1, -1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, -1, 1, 1, 0]
c40 = c20 + [-x for x in c20]
A20 = [2**(20 - (j + 2)//2) * c40[j] for j in range(40)]      # a_j = 2^(20 - ceil(j/2)) c_j, j = 1..40
rep('P(1,16,8,(8,8,4,0,-2,-2,-1,0)) - 4 pi', P(1, 16, 8, [8, 8, 4, 0, -2, -2, -1, 0], 60) - 4*pi)
rep('P(1,64,12,A8) - 8 sqrt2 pi', P(1, 64, 12, [32, 0, 8, 0, 8, 0, -4, 0, -1, 0, -1, 0], 40) - 8*sqrt(2)*pi)
rep('P(1,2^20,40,A) - 2^19 pi sqrt5/5', P(1, 2**20, 40, A20, 12) - 2**19*pi*sqrt(5)/5)
# T1, T2 of Section 5
C1, C2 = Cl2(2*pi/5), Cl2(4*pi/5)
z = lambda a: exp(mpc(0, 1)*pi*a/20) / sqrt(2)
ImL = lambda a: im(polylog(2, z(a)))
rep('T1 = ImLi2(z7)+ImLi2(z17) - [(21C1-C2)/16 - (3pi/20) log2]', ImL(7) + ImL(17) - ((21*C1 - C2)/16 - 3*pi/20*log(2)))
rep('T2 = ImLi2(z1)-ImLi2(z9) - [-(C1+21C2)/16 + (pi/20) log2]', ImL(1) - ImL(9) - (-(C1 + 21*C2)/16 + pi/20*log(2)))
# Proposition 5.1 at 600 digits
mp.dps = 620
S1 = (zeta(2, mpf(1)/5) - zeta(2, mpf(4)/5)) / 25
S2 = (zeta(2, mpf(2)/5) - zeta(2, mpf(3)/5)) / 25
al, be = sqrt(10 + 2*sqrt(5)), sqrt(10 - 2*sqrt(5))
V1 = [11010048, 0, -262144, -5505024, 0, 0, -65536, -65536, -688128, 0, -344064, -16384, -8192, 0, 0, -86016, -2048, 0,
      21504, 0, -10752, 0, 256, 5376, 0, 0, 64, 64, 672, 0, 336, 16, 8, 0, 0, 84, 2, 0, -21, 0]
V2 = [524288, 0, 5505024, -262144, 0, 0, 1376256, 1376256, -32768, 0, -16384, 344064, 172032, 0, 0, -4096, 43008, 0,
      1024, 0, -512, 0, -5376, 256, 0, 0, -1344, -1344, 32, 0, 16, -336, -168, 0, 0, 4, -42, 0, -1, 0]
N = 2**16 * 221
rep('Prop 5.1 (i):  S1 - (4 pi log2/5525)(32a+9b) - P(V1)/(2^16*221)', S1 - 4*pi*log(2)/5525*(32*al + 9*be) - P(2, 2**20, 40, V1, 110)/N)
rep('Prop 5.1 (ii): S2 + (4 pi log2/5525)(9a-32b) - P(V2)/(2^16*221)', S2 + 4*pi*log(2)/5525*(9*al - 32*be) - P(2, 2**20, 40, V2, 110)/N)
rep('red control: (i) with 32 replaced by 31', S1 - 4*pi*log(2)/5525*(31*al + 9*be) - P(2, 2**20, 40, V1, 110)/N)
print('\n'.join(out))
