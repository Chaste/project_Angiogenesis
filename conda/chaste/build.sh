
mkdir build
cd build

#CC=$PREFIX/bin/cc
#CXX=$PREFIX/bin/c++
export LIBRARY_PATH=$PREFIX/lib
export INCLUDE_PATH=$PREFIX/include

cmake .. \
    -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DBUILD_TESTING:BOOL=OFF \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DBOOST_ROOT=$PREFIX \
    -DChaste_ENABLE_TESTING=ON \
    -DChaste_UPDATE_PROVENANCE=OFF \
    -DChaste_ENABLE_heart_TESTING=OFF \
    -DChaste_ENABLE_lung_TESTING=OFF \
    -DChaste_ENABLE_crypt_TESTING=OFF \
    -DChaste_ENABLE_global_TESTING=OFF \
    -DChaste_ENABLE_linalg_TESTING=OFF \
    -DChaste_ENABLE_io_TESTING=OFF \
    -DChaste_ENABLE_mesh_TESTING=OFF \
    -DChaste_ENABLE_ode_TESTING=OFF \
    -DChaste_ENABLE_pde_TESTING=OFF \
    -DChaste_ENABLE_cell_based_TESTING=OFF \
    -DChaste_ENABLE_continuum_mechanics_TESTING=OFF \
    -DChaste_ENABLE_project_Angiogenesis_TESTING=OFF \
    -DVTK_DIR=$PREFIX \
    -DXERCESC_LIBRARY=$LIBRARY_PATH/libxerces-c.so \
    -DXERCESC_INCLUDE=$INCLUDE_PATH \
    -DXSD_EXECUTABLE=$PREFIX/bin/xsd
#    -DHDF5_DIR=${HDF5_DIR} 
#    -DMPI_DIR=${MPI_DIR} 
#    -DPETSC_DIR=${PETSC_DIR} 
#    -DPETSC_ARCH=${PETSC_ARCH}
#     -SPYTHON_INCLUDE_DIR:FILEPATH=$PYTHON \
make chaste_project_Angiogenesis -j $CPU_COUNT
make project_Angiogenesis_Python -j $CPU_COUNT
make install -j $CPU_COUNT
cp projects/Angiogenesis/libchaste_project_Angiogenesis.so $LIBRARY_PATH/chaste/
cd projects/Angiogenesis/python
python setup.py install --prefix=$PREFIX

