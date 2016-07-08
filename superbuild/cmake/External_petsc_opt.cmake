# Build petsc

list(APPEND projects petsc_OPT)
set(_v 6.3)
set(petsc_version 3.${_v})
set(petsc_url "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${petsc_version}.tar.gz")

# Need openmpi libraries to be in path, DOING MANUALLY FOR NOW
set(ENV{LD_LIBRARY_PATH} "${MPI_DIR}/lib:$ENV{LD_LIBRARY_PATH}")

set(petsc_cmds
CONFIGURE_COMMAND <SOURCE_DIR>/config/configure.py PETSC_DIR=<SOURCE_DIR> PETSC_ARCH=linux-gnu-opt --download-f2cblaslapack --with-mpi-dir=${MPI_DIR} --download-hypre --download-scalapack --download-ptscotch=1 --download-mumps --download-sundials --download-hdf5 --download-superlu_dist --download-parmetis --download-metis --with-x=false --with-clanguage=cxx --with-shared-libraries --with-debugging=0 --download-suitesparse --prefix=<INSTALL_DIR>
BUILD_COMMAND make PETSC_DIR=<SOURCE_DIR> PETSC_ARCH=linux-gnu-opt
)

ExternalProject_Add(petsc_opt
  DEPENDS openmpi
  URL ${petsc_url}
  ${petsc_cmds}
  INSTALL_COMMAND make install DESTDIR=<INSTALL_DIR>
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(petsc_opt source_dir)
set(HDF5_DIR "${source_dir}/linux-gnu-opt/" CACHE INTERNAL "")

ExternalProject_Get_Property(petsc_opt source_dir)
set(PETSC_DIR "${source_dir}" CACHE INTERNAL "")
set(PETSC_ARCH "linux-gnu-opt" CACHE INTERNAL "")
