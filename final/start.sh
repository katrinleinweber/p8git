
lammps=/home/va83ruvu/tallylammps/src/lmp_mpi
N=$(pwd|awk -F"/" '{print $NF}')
folder=$(pwd|awk -F"/" '{print $(NF-1)}')

a=4.08;
export OMP_NUM_THREADS=1

/home/va83ruvu/tallylammps/src/lmp_mpi -var a $a -var N $N -in lattice.in
/home/va83ruvu/tallylammps/src/lmp_mpi -var a $a -in create.in

K=10;
perl -pi -e s/PROC/4/g job.sh
perl -pi -e s/NAME/$folder.$N/g job.sh 

#qsub job.sh $a $K
