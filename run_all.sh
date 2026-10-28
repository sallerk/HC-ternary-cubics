#!/bin/sh
# run_all.sh -- rerun the computations of this repository.
# Usage:   sh run_all.sh [OUT]
#   OUT  the folder for the rerun (default: rerun). The scripts are copied to OUT/code and write their results
#        into OUT/outputs, so the reference outputs in outputs/ are not touched; their screen output goes to
#        OUT/logs/NAME.log.
# Then:    python3 -I compare_outputs.py outputs OUT/outputs
# The interpreter is python3, or the value of the environment variable PYTHON (for example PYTHON=python).
# The scripts run one after the other, in this order, because t3 reads the output t2_hodge_count.json of t2
# (see README.txt, Sections 2 and 5). A script that fails does not stop the others; each line "NAME exit N"
# reports the exit status of one run.
OUT=${1:-rerun}
PY=${PYTHON:-python3}
T=$(dirname "$0")
C="$OUT/code"
O="$OUT/outputs"
L="$OUT/logs"
mkdir -p "$C" "$O" "$L" || exit 1
cp "$T"/code/*.py "$C"/ || exit 1
run() {
  name=$1
  "$PY" -I -B "$C/$name.py" > "$L/$name.log" 2>&1
  echo "$name exit $?"
}
run t1_sk_smooth
run t2_hodge_count
run t3_formula
run t4_tenfold_range
run t5_quaternary
run t6_joins
run t7_eckardt
run t7_quinary
