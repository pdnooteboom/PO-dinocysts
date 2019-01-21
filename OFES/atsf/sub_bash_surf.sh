#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -N 1

for i in {0..49}
do
python atsf_OFES_global_noadv.py $i &
sleep 5
done

# wait until all background processes are ended:

wait
