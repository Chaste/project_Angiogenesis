# Build UFL

list(APPEND projects fiat)

ExternalProject_Add(fiat
  GIT_REPOSITORY https://bitbucket.org/fenics-project/fiat.git
  GIT_TAG "fiat-1.6.0"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND <SOURCE_DIR>/setup.py install --prefix=/scratch/jgrogan/Software/miniconda2/
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(fiat install_dir)
set(FIAT_DIR "${install_dir}" CACHE INTERNAL "")
