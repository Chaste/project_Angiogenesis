# Build dolfin

list(APPEND projects dolfin)

ExternalProject_Add(dolfin
  DEPENDS boost petsc ffc vtk
  GIT_REPOSITORY https://bitbucket.org/fenics-project/dolfin.git
  GIT_TAG ""
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> 
        -DBOOST_ROOT=${BOOST_ROOT}
        -DBOOST_INCLUDEDIR=${BOOST_ROOT}/include
        -DPETSC_DIR=${CMAKE_CURRENT_SOURCE_DIR}/external/Install/petsc_opt
        -DPETSC_ARCH=
        -DVTK_DIR=${VTK_DIR} 
        -DMPIEXEC=${MPI_DIR}/bin/mpiexec
        -DMPI_DIR=${MPI_DIR}
        -DMPI_INCLUDE_DIR=${MPI_DIR}/include
        -DMPI_DIR=${MPI_DIR}
        -DMPI_C_COMPILER=${MPI_DIR}/bin/mpicc
        -DMPI_CXX_COMPILER=${MPI_DIR}/bin/mpicxx
        -DDOLFIN_ENABLE_SPHINX=0FF
        -DDOLFIN_ENABLE_SLEPC=0FF
        -DDOLFIN_ENABLE_TILINOS=0FF
        -DDOLFIN_ENABLE_HDF5=0FF
        -DDOLFIN_ENABLE_HDF5=0FF
        -DUFC_DIR=/home/grogan/.pyenv/versions/2.7.11/share/ufc  
        -DMPI_CXX_COMPILER=${MPI_DIR}/bin/mpicxx
    BUILD_COMMAND make
    INSTALL_COMMAND make install
)
