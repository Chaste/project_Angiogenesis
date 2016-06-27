# Build UFL

list(APPEND projects ufl)

ExternalProject_Add(ufl
  GIT_REPOSITORY https://bitbucket.org/fenics-project/ufl.git
  GIT_TAG "ufl-1.6.0"
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND <SOURCE_DIR>/setup.py install --prefix=/scratch/jgrogan/Software/miniconda2/
  BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(ufl install_dir)
set(UFL "${install_dir}" CACHE INTERNAL "")
