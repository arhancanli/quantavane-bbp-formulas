#!/bin/zsh
# Sweep 7 (11 Nov 2026 23:1x, coordinator): (A) base 5, weight 2, the Q(sqrt5)/Q(sqrt-5) family with L(2,chi_-20) and
# L(2,chi_-15); (B) base 2^3 with period 7, weight 2, the Q(sqrt-7) family with L(2,chi_-7) (the next field after our
# proved Q(sqrt-2) formula). Any hit is only a candidate: it is re-checked at higher precision and then must be PROVED.
cd "${0:A:h}/.."
TA='pi**2;l5**2;lg**2;l5*lg;s5*pi**2;s5*lg**2;s5*l5*lg;L(-20,2);s5*L(-20,2);L(-15,2);s15*L(-15,2);s3*L(-15,2);G;pi*l5;pi*lg;s5*pi*lg;s5*pi*l5'
echo "### (A) base 5, weight 2, Q(sqrt5) family"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 5 2 "1:4,2:8,1:10,2:10,1:20,2:20,4:20" "$TA"
TB='pi**2;l2**2;pi*l2;L(-7,2);s7*L(-7,2);pi*l7;s7*pi*l7;l7**2;l2*l7;s7*pi**2;s7*pi*l2;G;L(-8,2);L(-3,2);s3*L(-3,2)'
echo "### (B) base 2^3, period 7, weight 2, Q(sqrt-7) family"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 2 2 "3:7,3:14,6:14,3:21,6:21,3:28,6:28,6:42" "$TB"
echo "SWEEP7 ALL DONE"
