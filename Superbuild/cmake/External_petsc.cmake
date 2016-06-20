# Build petsc

list(APPEND projects petsc)
set(_v 6.3)
set(petsc_version 3.${_v})
set(petsc_url "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-${petsc_version}.tar.gz")

# Need openmpi libraries to be in path, DOING MANUALLY FOR NOW
set(ENV{LD_LIBRARY_PATH} "${MPI_DIR}/lib:$ENV{LD_LIBRARY_PATH}")

set(petsc_cmds
CONFIGURE_COMMAND <SOURCE_DIR>/config/configure.py PETSC_DIR=<SOURCE_DIR> PETSC_ARCH=linux-gnu --download-f2cblaslapack --with-mpi-dir=${MPI_DIR} --download-hypre --download-sundials --download-hdf5 --download-parmetis --download-metis --with-x=false --with-clanguage=cxx --with-shared-libraries
BUILD_COMMAND make PETSC_DIR=<SOURCE_DIR> PETSC_ARCH=linux-gnu
)

ExternalProject_Add(petsc
  DEPENDS openmpi
  URL ${petsc_url}
  ${petsc_cmds}
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(petsc source_dir)
set(HDF5_DIR "${source_dir}/linux-gnu/" CACHE INTERNAL "")

ExternalProject_Get_Property(petsc source_dir)
set(PETSC_DIR "${source_dir}" CACHE INTERNAL "")
set(PETSC_ARCH "linux-gnu" CACHE INTERNAL "")
