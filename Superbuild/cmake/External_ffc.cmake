# Build UFL

list(APPEND projects ffc)

set(ENV{PATH} "${SWIG_DIR}/bin:$ENV{PATH}")

ExternalProject_Add(ffc
  DEPENDS ufl fiat instant swig
  GIT_REPOSITORY https://bitbucket.org/fenics-project/ffc.git
  GIT_TAG ""
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND export "PATH=${SWIG_DIR}/bin:$PATH"
  BUILD_COMMAND python <SOURCE_DIR>/setup.py install
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
)
