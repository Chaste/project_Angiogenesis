# Build PETSC4PY

list(APPEND projects petsc4py)

set(ENV{PETSC_DIR} "${PETSC_DIR}")
set(ENV{PETSC_ARCH} "${PETSC_ARCH}")

ExternalProject_Add(petsc4py
  DEPENDS petsc
  GIT_REPOSITORY https://bitbucket.org/petsc/petsc4py.git
  GIT_TAG "3.6.0" # careful, must match petsc version, e.g. 3.6
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND python <SOURCE_DIR>/setup.py install
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(petsc4py install_dir)
set(PETSC4PY_DIR "${install_dir}" CACHE INTERNAL "")
