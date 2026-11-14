"""Sweep 4: degree-2 BBP spans in mixed bases (6, 10, 12, 18, 24, 5, 7) -- any formula at all?"""
import sys, time
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from sweep3 import run
T = ['pi**2', 'l2**2', 'l3**2', 'l2*l3', 'pi*l2', 'pi*l3', 's3*pi*l2', 's3*pi*l3', 'G', 'L(-3,2)',
     'l5**2', 'l2*l5', 'l3*l5', 'pi*l5', 'l7*l2', 'l7*l3', 'pi*l7', 's5*pi*l5', 's2*pi*l2', 's6*pi*l2', 's6*pi*l3']
for (p, m, n) in [(6,1,2),(6,2,4),(6,6,12),(6,3,6),(6,12,24),(6,6,24),(12,1,2),(12,2,4),(12,6,12),(12,3,12),
                  (10,1,2),(10,2,4),(10,5,10),(10,10,20),(10,5,20),(5,2,4),(5,5,10),(5,10,20),(7,3,6),(7,6,12),
                  (18,2,4),(18,6,12),(24,2,4),(24,6,12)]:
    run(p, 2, m, n, T)
print("SWEEP4 DONE", flush=True)
