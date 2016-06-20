# Build UFL

list(APPEND projects ufl)

ExternalProject_Add(ufl
  GIT_REPOSITORY https://bitbucket.org/fenics-project/ufl.git
  GIT_TAG ""
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND python <SOURCE_DIR>/setup.py install
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
)
