#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -N 1 

python totimeseries_per_location_latlonadv_netcdf_cop.py

# wait until all background processes are ended:

wait
