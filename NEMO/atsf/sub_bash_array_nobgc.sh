#!/bin/bash
#SBATCH -t 2-02:00:00
##SBATCH -N 6 --ntasks-per-node=9
##SBATCH --mem-per-cpu=10000
#SBATCH --array=0-8 # starts 11 array tasks

# Job offset: <array_id> * 3 to <array_id>*3 + 2
START=$((SLURM_ARRAY_TASK_ID * 2))
END=$((START + 1))

# Run from START to END
for i in $(seq $START $END)
do
    python atsf_cart_nobgc.py  $i > out/atsf_cart_$i.out &
    sleep 5
done


wait
