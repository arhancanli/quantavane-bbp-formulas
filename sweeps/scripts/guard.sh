#!/bin/zsh
# Wait while the 1-minute load average is above 48 (formulas lane rule).
while true; do
  l=$(sysctl -n vm.loadavg | awk '{print $2}')
  if (( $(echo "$l < 48" | bc -l) )); then break; fi
  echo "[guard] load $l >= 48, waiting 60 s" >&2
  sleep 60
done
