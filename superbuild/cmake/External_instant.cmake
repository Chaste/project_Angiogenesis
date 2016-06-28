# Build UFL

list(APPEND projects instant)

ExternalProject_Add(instant
  GIT_REPOSITORY https://bitbucket.org/fenics-project/instant.git
  GIT_TAG "instant-1.6.0"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND <SOURCE_DIR>/setup.py install --prefix=/scratch/jgrogan/Software/miniconda2/
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(instant install_dir)
set(INSTANT_DIR "${install_dir}" CACHE INTERNAL "")
