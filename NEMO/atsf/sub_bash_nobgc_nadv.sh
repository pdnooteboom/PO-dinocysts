#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH -N 1
#SBATCH --mem-per-cpu=20000

for i in {0..17}
do
python atsf_cart_nobgc_nadv.py $i &
sleep 5
done

# wait until all background processes are ended:

wait
