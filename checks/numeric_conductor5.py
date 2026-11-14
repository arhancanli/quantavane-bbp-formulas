# Numerical check (not a proof) of the two conductor-5 binary formulas reported as observations in the paper.
from mpmath import mp, mpf, sqrt, pi, log, zeta, nstr
mp.dps = 600
S1 = (zeta(2, mpf(1)/5) - zeta(2, mpf(4)/5)) / 25          # sum 1/(5n+1)^2 - 1/(5n+4)^2
S2 = (zeta(2, mpf(2)/5) - zeta(2, mpf(3)/5)) / 25          # sum 1/(5n+2)^2 - 1/(5n+3)^2
al, be = sqrt(10 + 2*sqrt(5)), sqrt(10 - 2*sqrt(5))
A1 = [13893632, -786432, -393216, -6029312, 786432, -196608, -65536, -491520, -851968, -49152, -434176, 32768, -12288, -12288, -24576, -133120, -2048, -3072, 26624, -3072, -13568, -768, 384, 5120, -768, -192, 64, 288, 832, -48, 424, -80, 12, -12, 24, 118, 2, -3, -26, 0]
A2 = [786432, -786432, 6815744, 524288, 786432, -196608, 1736704, 1310720, -32768, -49152, -24576, 483328, 212992, -12288, -24576, -30720, 54272, -3072, 1024, -3072, -768, -768, -6656, -1280, -768, -192, -1696, -1472, 32, -48, 24, -520, -208, -12, 24, 18, -53, -3, -1, 0]
P = lambda A: sum(mpf(2)**(-20*k) * sum(mpf(A[j]) / (40*k + j + 1)**2 for j in range(40)) for k in range(0, 110))
print('S1 =', nstr(S1, 25), ' S2 =', nstr(S2, 25))
print('|5 S1 - (pi log2/1105)(128 al + 36 be) - P(A1)/3620864| =', nstr(5*S1 - pi*log(2)/1105*(128*al + 36*be) - P(A1)/3620864, 5))
print('|5 S2 + (pi log2/1105)(36 al - 128 be) - P(A2)/3620864| =', nstr(5*S2 + pi*log(2)/1105*(36*al - 128*be) - P(A2)/3620864, 5))
from mpmath import clsin
print('|Cl2(2pi/5) - (al S1 + be S2)/4| =', nstr(clsin(2, 2*pi/5) - (al*S1 + be*S2)/4, 5), ' |Cl2(4pi/5) - (be S1 - al S2)/4| =', nstr(clsin(2, 4*pi/5) - (be*S1 - al*S2)/4, 5))
