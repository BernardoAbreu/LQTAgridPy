#!/bin/bash

# declare -a array=("81Mols" "antimalaricos" "autismo" "cance_prostata" "hiv" "TP")
declare -a array=("TP")

for DATA in "${array[@]}"; do
	# for THETA in 2 3 4 6 8 10; do
	for THETA in 10; do
		mkdir $DATA
		echo "Running ${DATA} dataset with theta=${THETA}"
		python -u src/lqtagrid_2.py --mols "Datasets/${DATA}" -a NH3+ -r 1.0 -d $THETA -i 2.5 -o "${DATA}/d${THETA}"
	done
done
