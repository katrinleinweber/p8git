#!/bin/bash
#SBATCH --mem-per-cpu=100
#SBATCH -n 4
#SBATCH -C avx2
#SBATCH -t 9:59:59
#SBATCH --job-name="size.20"
 
n=$(ls|grep run|wc -l)
n=$(expr $n + 1)
echo "hey">$n.run
#F=$(pwd| awk -F"/" '{print $(NF-1)}');
eps=$(pwd| awk -F"/" '{print $(NF)}');


srun /home/va83ruvu/tallylammps/src/lmp_mpi -var n $n -var a $1 -var K $2 -in test.in>test_$n.lammps
