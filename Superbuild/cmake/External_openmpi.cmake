# Build openmpi

list(APPEND projects openmpi)
set(openmpi_url "https://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.8.tar.bz2")

set(openmpi_cmds
CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
BUILD_COMMAND make all
INSTALL_COMMAND make install
)

ExternalProject_Add(openmpi
  URL ${openmpi_url}
  ${openmpi_cmds}
)

ExternalProject_Get_Property(openmpi install_dir)
set(MPI_DIR "${install_dir}" CACHE INTERNAL "")
