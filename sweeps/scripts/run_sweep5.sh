#!/bin/zsh
# Sweep 5: chi_-24 in binary; chi_-8 and chi_-24 in ternary (degree 2), with natural partners.
cd "${0:A:h}/.."
TB='pi**2;G;pi*l2;l2**2;L(-3,2);L(-8,2);s2*pi*l2;s3*pi*l2;L(-24,2);s6*L(-24,2);s2*L(-24,2);s3*L(-24,2);s6*pi*l2;s6*pi*l3;s2*pi*l3;pi*l3;s3*pi*l3;l2*l3;l3**2;s6*pi*lt;pi*lt;s2*pi*ls;s3*pi*ls'
echo "### binary chi_-24"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 2 2 "12:24,6:24,24:48,12:48" "$TB"
TT='pi**2;l3**2;s3*pi*l3;L(-3,2);G;pi*l3;pi*l2;L(-8,2);s2*L(-8,2);s2*pi*l3;s2*pi*l2;s6*pi*l3;s6*pi*l2;L(-24,2);s6*L(-24,2);s2*L(-24,2);s3*L(-24,2);l2*l3;l2**2;s3*pi*l2;s2*pi*ls;s6*pi*lt;pi*lt'
echo "### ternary chi_-8, chi_-24"
./scripts/guard.sh; nice -n 10 ./venv/bin/python -I scripts/sweep3.py 3 2 "12:24,6:24,24:24,24:48,12:48" "$TT"
echo "SWEEP5 ALL DONE"
