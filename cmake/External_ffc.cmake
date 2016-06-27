# Build UFL

list(APPEND projects ffc)

set(ENV{PATH} "${SWIG_DIR}/bin:$ENV{PATH}")

ExternalProject_Add(ffc
  DEPENDS ufl fiat instant swig
  GIT_REPOSITORY https://bitbucket.org/fenics-project/ffc.git
  GIT_TAG "ffc-1.6.0"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND export "PATH=${SWIG_DIR}/bin:$PATH"
  BUILD_COMMAND ""
  INSTALL_COMMAND python <SOURCE_DIR>/setup.py install --prefix=/scratch/jgrogan/Software/miniconda2/
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(ffc install_dir)
set(FFC_DIR "${install_dir}" CACHE INTERNAL "")
set(UFC_DIR "${install_dir}" CACHE INTERNAL "")
