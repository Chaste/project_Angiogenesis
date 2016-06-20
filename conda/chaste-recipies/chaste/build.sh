
mkdir build
cd build

CC=$PREFIX/bin/cc
CXX=$PREFIX/bin/c++

export LIBRARY_PATH=$PREFIX/lib
export INCLUDE_PATH=$PREFIX/include

export BLAS_DIR=$LIBRARY_PATH

cmake .. \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DBUILD_TESTING:BOOL=OFF \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DBOOST_ROOT=$PREFIX \
    -DChaste_ENABLE_TESTING=ON \
    -DChaste_UPDATE_PROVENANCE=OFF \
    -DChaste_ENABLE_heart_TESTING=OFF \
    -DChaste_ENABLE_lung_TESTING=OFF \
    -DChaste_ENABLE_crypt_TESTING=OFF \
    -DChaste_ENABLE_continuum_mechanics_TESTING=OFF

#    -DHDF5_DIR=${HDF5_DIR} 
#    -DMPI_DIR=${MPI_DIR} 
#    -DPETSC_DIR=${PETSC_DIR} 
#    -DPETSC_ARCH=${PETSC_ARCH}
#    -DVTK_DIR=$PREFIX \
#    -DXERCES_LIBRARY=${LIBRARY_PATH} \
#    -DXERCES_INCLUDE=${INCLUDE_PATH} \
#    -DBOOST_ROOT=$PREFIX \
#    -DXSD_EXECUTABLE=$PREFIX/bin/xsd \
#    -DXSD_INCLUDE_DIR=$PREFIX/libxsd/xsd/cxx \

make chaste_project_Angiogenesis -j 3
make project_Angiogenesis_Python
python projects/Angiogenesis/python/setup.py install
#make install
