#!/bin/sh
#SBATCH --time=00:30:00
#SBATCH -N 1 -n 4
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --job-name=weibel

# specify all the variables

EXECUTABLE=tristan-mp2d
INPUT=input.weibel
OUTPUT_DIR=output2
SLICE_DIR=slices
RESTART_DIR=restart
REPORT_FILE=report
ERROR_FILE=error


module load intel-rt/2021.1.2
module load intel-tbb/2021.1.1
module load intel-mkl/2021.1.1
module load intel-debugger/10.0.0
module load intel-dpl/2021.1.2
module load /opt/intel/oneapi/compiler/2021.1.2/linux/lib/oclfpga/modulefiles/oclfpga
module load intel/2021.1.2
module load ucx/1.9.0
module load intel-mpi/intel/2021.1.1
module load hdf5/intel-2021.1/intel-mpi/1.10.6
module load anaconda3/2020.11

# create the output directory

mkdir $OUTPUT_DIR
# backup the executable and the input file

cp $EXECUTABLE $OUTPUT_DIR
cp $INPUT $OUTPUT_DIR
srun -N 1 -n 4 --ntasks-per-node=96 $EXECUTABLE -i $INPUT -o $OUTPUT_DIR -s $SLICE_DIR -r $RESTART_DIR > $OUTPUT_DIR/$REPORT_FILE 2> $OUTPUT_DIR/$ERROR_FILE
