## The GCMC-vasp code was run on High Performance Computing platforms with slurm.

## Slurm commands:

---------------------------------------------------------------------------------

cd $SLURM_SUBMIT_DIR

echo "Job started on `hostname` at `date`"

module purge
module load intel/2018b
module load ifort/2018.4.274-GCC-7.3.0  impi/2018.4.274 netCDF-Fortran/4.4.4

ifort vasp_FeC_system.f90 -o  GCMC_vasp ;

./GCMC_vasp

rm vasp_std

echo "Job finished on `hostname` at `date`"

---------------------------------------------------------------------------------
