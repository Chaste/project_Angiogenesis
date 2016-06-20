# Build chaste

list(APPEND projects chaste)

set(chaste_cmds
CMAKE_ARGS -DBUILD_SHARED_LIBS:BOOL=ON 
    -DBUILD_TESTING:BOOL=OFF 
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> 
    -DChaste_ENABLE_TESTING=ON 
    -DHDF5_DIR=${HDF5_DIR} 
    -DMPI_DIR=${MPI_DIR} 
    -DPETSC_DIR=${PETSC_DIR} 
    -DPETSC_ARCH=${PETSC_ARCH}
    -DVTK_DIR=${VTK_DIR} 
    -DXERCES_LIBRARY=${XERCES_DIR}/lib 
    -DXERCES_INCLUDE=${XERCES_DIR}/include 
    -DBOOST_ROOT=${BOOST_ROOT} 
    -DXSD_EXECUTABLE=${XSD_DIR}/bin/xsd
    -DXSD_INCLUDE_DIR=${XSD_DIR}/libxsd/xsd/cxx
    -DChaste_UPDATE_PROVENANCE=OFF
    -DChaste_ENABLE_heart_TESTING=OFF  
    -DChaste_ENABLE_lung_TESTING=OFF
    -DChaste_ENABLE_crypt_TESTING=OFF
    -DChaste_ENABLE_continuum_mechanics_TESTING=OFF

PATCH_COMMAND patch -t -N < ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Patches/patch_chaste_make.patch
BUILD_COMMAND make chaste_cell_based
)

ExternalProject_Add(chaste
  DEPENDS boost petsc xsd xerces vtk
  SVN_REPOSITORY https://chaste.cs.ox.ac.uk/svn/chaste/trunk/ 
  ${chaste_cmds}
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(chaste install_dir)
set(CHASTE_DIR "${install_dir}" CACHE INTERNAL "")
