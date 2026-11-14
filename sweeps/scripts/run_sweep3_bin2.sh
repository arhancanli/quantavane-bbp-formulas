#!/bin/zsh
cd "${0:A:h}/.."
T='pi**2;G;pi*l2;l2**2;s3*L(-3,2);L(-3,2);L(-8,2);s2*L(-8,2);L(-7,2);s7*L(-7,2);s2*pi*l2;s3*pi*l2;s2*pi*ls;pi*ls;l2*ls;s2*l2*ls;ls**2;l3**2;l2*l3;pi*l3;s3*pi*l3;s2*pi*l3;L(-20,2);s5*L(-20,2);L(-24,2);s6*L(-24,2)'
./scripts/guard.sh
nice -n 10 ./venv/bin/python -I scripts/sweep3.py 2 2 "12:24,24:24,6:24,24:48,12:48,48:48,20:40,10:40,40:40,16:32,8:32,32:32" "$T"
