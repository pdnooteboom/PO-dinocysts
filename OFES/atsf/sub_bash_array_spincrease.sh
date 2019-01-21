#!/bin/bash
#SBATCH -t 4-00:00:00
#SBATCH -N 1 --ntasks-per-node=5
##SBATCH --mem-per-cpu=20000
#SBATCH --array=0-4 # starts 11 array tasks

echo SLURM_ARRAY_TASK_ID
# Job offset: <array_id> * 3 to <array_id>*3 + 2
START=$((SLURM_ARRAY_TASK_ID * 10))
END=$((START + 9))

# Run from START to END
for i in $(seq $START $END)
do
    python atsf_OFES_global_spincrease_depth.py $i &
    sleep 5
done


wait
