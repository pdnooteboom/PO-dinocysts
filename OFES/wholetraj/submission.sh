#/bin/sh
# SGE: the job name
#$ -N one_particle_adv
#
# The requested run-time, expressed as (xxxx sec or hh:mm:ss)
#$ -l h_rt=40:00:00
#
# SGE: your Email here, for job notification
#$ -M p.d.nooteboom@uu.nl
#
# SGE: when do you want to be notified (b : begin, e : end, s : error)?
#$ -m e
#
# SGE: ouput in the current working dir
#$ -cwd    
#

cd /home/students/3830241/PhD/OFES_global/wholetraj/
python loc_traj_OFES.py -p 10 -v
