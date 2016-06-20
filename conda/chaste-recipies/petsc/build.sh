#!/bin/bash

export LIBRARY_PATH=$PREFIX/lib

./configure \
  --prefix=$PREFIX \
  --with-mpi-dir=$PREFIX \
  --download-suitesparse \
  --download-hypre \
  --download-f2cblaslapack \
  --download-hdf5 \
  --download-parmetis \
  --download-metis \
  --download-sundials \
  --with-x=false \
  --with-clanguage=cxx \
  --with-shared-libraries
make
make install

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
